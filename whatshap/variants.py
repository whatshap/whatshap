"""
Detect variants in reads.
"""
import logging
from collections import defaultdict, Counter, deque
from typing import Iterable, Iterator, List, Optional
from dataclasses import dataclass

from .core import Read, ReadSet, NumericSampleIds
from .bam import SampleBamReader, MultiBamReader, BamReader
from .align import edit_distance, edit_distance_affine_gap
from ._variants import _iterate_cigar


logger = logging.getLogger(__name__)


class ReadSetError(Exception):
    pass


@dataclass
class AlleleProgress:
    progress: int = 0
    length: int = 0
    quality: int = 0
    matched: int = 0
    match_target: int = 0
    inserted: int = 0
    insert_target: int = 0
    deleted: int = 0
    delete_target: int = 0


class VariantProgress:
    def __init__(self, variant_id):
        self.variant_id = variant_id
        self.query_start = 0
        self.alleles = []

    def __iter__(self):
        for a in self.alleles:
            yield a

    def __len__(self):
        return len(self.alleles)

    def add_allele(self, matches, insertions, deletions):
        l = matches + insertions + deletions
        a = AlleleProgress(0, l, 0, 0, matches, 0, insertions, 0, deletions)
        self.alleles.append(a)

    def reset(self, query_start):
        self.query_start = query_start
        for a in self.alleles:
            a.progress, a.matched, a.inserted, a.deleted, a.quality = 0, 0, 0, 0, 0

    def get_resolved(self):
        return [i for i, a in enumerate(self.alleles) if a.progress == a.length]

    def get_pending(self):
        return [i for i, a in enumerate(self.alleles) if 0 <= a.progress < a.length]


class ReadSetReader:
    """
    Associate VCF variants with BAM reads.

    A VCF file contains variants, and a BAM file contain reads, but the
    information which read contains which variant is not available. This
    class re-discovers the variants in each read, using the
    knowledge in the VCF of where they should occur.
    """

    def __init__(
        self,
        paths: List[str],
        reference: Optional[str],
        numeric_sample_ids: NumericSampleIds,
        mapq_threshold: int = 20,
        overhang: int = 10,
        affine: int = False,
        gap_start: int = 10,
        gap_extend: int = 7,
        default_mismatch: int = 15,
    ):
        """
        paths -- list of BAM paths
        reference -- path to reference FASTA (can be None)
        numeric_sample_ids -- ??
        mapq_threshold -- minimum mapping quality
        overhang -- extend alignment by this many bases to left and right
        affine -- use affine gap costs
        gap_start, gap_extend, default_mismatch -- parameters for affine gap cost alignment
        """
        self._mapq_threshold = mapq_threshold
        self._numeric_sample_ids = numeric_sample_ids
        self._use_affine = affine
        self._gap_start = gap_start
        self._gap_extend = gap_extend
        self._default_mismatch = default_mismatch
        self._overhang = overhang
        self._paths = paths
        self._reader: BamReader
        if len(paths) == 1:
            self._reader = SampleBamReader(paths[0], reference=reference)
        else:
            self._reader = MultiBamReader(paths, reference=reference)

    @property
    def n_paths(self):
        return len(self._paths)

    def read(self, chromosome, variants, sample, reference, regions=None) -> ReadSet:
        """
        Detect alleles and return a ReadSet object containing reads representing
        the given variants.

        If a reference is provided (reference is not None), alleles are
        detected by re-aligning sections of the query to the REF and ALT
        sequence extended a few bases to the left and right.

        If reference is None, alleles are detected by inspecting the
        existing alignment (via the CIGAR).

        chromosome -- name of chromosome to work on
        variants -- list of vcf.VcfVariant objects
        sample -- name of sample to work on. If None, read group information is
            ignored and all reads in the file are used.
        reference -- reference sequence of the given chromosome (or None)
        regions -- list of start,end tuples (end can be None)
        """
        # Since variants are identified by position, positions must be unique.
        if __debug__ and variants:
            varposc = Counter(variant.position for variant in variants)
            pos, count = varposc.most_common()[0]
            assert count == 1, f"Position {pos} occurs more than once in variant list."

        alignments = self._usable_alignments(chromosome, sample, regions)
        reads = self._alignments_to_reads(alignments, variants, sample, reference)
        grouped_reads = self._group_paired_reads(reads)
        readset = self._make_readset_from_grouped_reads(grouped_reads)
        return readset

    @staticmethod
    def _make_readset_from_grouped_reads(groups: Iterable[List[Read]]) -> ReadSet:
        read_set = ReadSet()
        for group in groups:
            read_set.add(merge_reads(*group))
        return read_set

    @staticmethod
    def _group_paired_reads(reads: Iterable[Read]) -> Iterator[List[Read]]:
        """
        Group reads into paired-end read pairs. Uses name, source_id and sample_id
        as grouping key.

        TODO
        Grouping by name should be sufficient since the SAM spec states:
        "Reads/segments having identical QNAME are regarded to come from the same template."
        """
        groups = defaultdict(list)
        for read in reads:
            groups[(read.source_id, read.name, read.sample_id)].append(read)
        for group in groups.values():
            if len(group) > 2:
                raise ReadSetError(
                    f"Read name {group[0].name!r} occurs more than twice in the input file"
                )
            yield group

    def _usable_alignments(self, chromosome, sample, regions=None):
        """
        Retrieve usable (suficient mapping quality, not secondary etc.)
        alignments from the alignment file
        """
        if regions is None:
            regions = [(0, None)]
        for s, e in regions:
            for alignment in self._reader.fetch(
                reference=chromosome, sample=sample, start=s, end=e
            ):
                # TODO handle additional alignments correctly!
                # find out why they are sometimes overlapping/redundant
                if (
                    alignment.bam_alignment.flag & 2048 != 0
                    or alignment.bam_alignment.mapping_quality < self._mapq_threshold
                    or alignment.bam_alignment.is_secondary
                    or alignment.bam_alignment.is_unmapped
                    or alignment.bam_alignment.is_duplicate
                ):
                    continue
                yield alignment

    def has_reference(self, chromosome):
        return self._reader.has_reference(chromosome)

    def _alignments_to_reads(self, alignments, variants, sample, reference):
        """
        Convert BAM alignments to Read objects.

        If reference is not None, alleles are detected through re-alignment.

        Yield Read objects.
        """
        # FIXME hard-coded zero
        numeric_sample_id = 0 if sample is None else self._numeric_sample_ids[sample]
        if reference is not None:
            # Copy the pyfaidx.FastaRecord into a str for faster access
            reference = reference[:]
            normalized_variants = variants
        else:
            normalized_variants = [variant.normalized() for variant in variants]

        i = 0  # index into variants

        # Mark overlapping variants before to make later variant handling more modular
        conflict_vars = self.detect_overlapping_variants(normalized_variants)

        # Create allele progress trackers once, instead of doing it for every read again
        if reference is None:
            var_progress = {
                j: self.build_var_progress(normalized_variants, j)
                for j in range(i, len(normalized_variants))  # if j not in conflict_vars
            }

        for alignment in alignments:
            # Skip variants that are to the left of this read
            while (
                i < len(normalized_variants)
                and normalized_variants[i].position < alignment.bam_alignment.reference_start
            ):
                i += 1

            try:
                barcode = alignment.bam_alignment.get_tag("BX")
            except KeyError:
                barcode = ""

            read = Read(
                alignment.bam_alignment.qname,
                alignment.bam_alignment.mapq,
                alignment.source_id,
                numeric_sample_id,
                alignment.bam_alignment.reference_start,
                barcode,
            )

            if reference is None:
                detected = self.detect_alleles(
                    normalized_variants, conflict_vars, var_progress, i, alignment.bam_alignment
                )
            else:
                detected = self.detect_alleles_by_alignment(
                    variants,
                    i,
                    alignment.bam_alignment,
                    reference,
                    self._overhang,
                    self._use_affine,
                    self._gap_start,
                    self._gap_extend,
                    self._default_mismatch,
                )

            for j, allele, quality in detected:
                read.add_variant(variants[j].position, allele, quality)
            if read:  # At least one variant covered and detected
                yield read

    def detect_alleles(self, variants, conflict_vars, var_progress, j, bam_read):
        """
        Detect the correct alleles of the variants that are covered by the
        given bam_read.

        Yield tuples (index, allele, quality), where index is into the variants list.

        variants -- list of variants (VcfVariant objects)
        j -- index of the first variant (in the variants list) to check
        """
        ref_pos = bam_read.reference_start  # position relative to reference
        query_pos = 0  # position relative to read

        # Skip variants that come before this region
        while j < len(variants) and variants[j].position < ref_pos:
            j += 1

        vqueue = deque()  # buffer for pending variants to keep them in positional order
        seen = set()  # seen variant positions

        for cigar_op, length in bam_read.cigartuples:
            # Skip variants that come before this region
            while j < len(variants) and (variants[j].position < ref_pos or j in conflict_vars):
                j += 1

            # MIDNSHPX= => 012345678. Skip for soft clipping/padding, etc.
            if cigar_op == 3:  # N operator (reference skip)
                ref_pos += length
                continue
            elif cigar_op == 4:  # S operator (soft clipping)
                query_pos += length
                continue
            elif cigar_op == 5 or cigar_op == 6:  # H or P (hard clipping or padding)
                continue

            # Queue all variants that start within the ref span of the cigar operation
            ref_end = ref_pos + length
            while j < len(variants) and variants[j].position < ref_end:
                # Skip overlapped variants
                if j in conflict_vars:
                    j += 1
                    continue

                # Skip duplicate positions (TODO: At least for now)
                if variants[j].position in seen:
                    j += 1
                    assert False
                    continue

                ref_len = len(variants[j].reference_allele)
                # Special case: If a non-insertion variant is seen by I-Op, continue with next Op
                # It is an insertion in front of a non-insertion variant, that must be ignored
                # We cannot skip I-Op in general, because this might overlook insertion variants
                if cigar_op == 1 and ref_len > 0:
                    break
                # Special case: If a D-Op sees an insertion variant, skip this variant to be conform
                # with old implementation. Actually, it would be correct to assume ref allele here,
                # if the preivous base matched. This seems to be an artifact of normalized variants.
                seen.add(variants[j].position)
                if cigar_op == 2 and ref_len == 0:
                    j += 1
                    continue

                # Entry = [var_id, query_start_pos, for each allele: AlleleProgress object]
                query_start = (
                    query_pos + variants[j].position - ref_pos if cigar_op != 2 else query_pos
                )
                var_progress[j].reset(query_start)
                vqueue.append(var_progress[j])
                j += 1

            # Handle detection and positional progress depending on cigar op
            ref_end = ref_pos
            query_end = query_pos
            if cigar_op in (0, 7, 8):  # M, X, = operators (match)
                handler = self.detect_alleles_match
                ref_end += length
                query_end += length
            elif cigar_op == 1:  # I operator (insertion)
                handler = self.detect_alleles_insertion
                query_end += length
            elif cigar_op == 2:  # D operator (deletion)
                handler = self.detect_alleles_deletion
                ref_end += length
            else:
                logger.error("Unsupported CIGAR operation: %d", cigar_op)
                raise ValueError(f"Unsupported CIGAR operation: {cigar_op}")

            # Progress pending variants using current cigar operation
            for var_entry in vqueue:
                variant = variants[var_entry.variant_id]
                handler(variant, var_entry, bam_read, ref_pos, query_pos, length)
            ref_pos = ref_end
            query_pos = query_end

            # Yield resolved variants from left, pop inresolvable variants
            while vqueue:
                var_entry = vqueue.popleft()
                resolved = list(var_entry.get_resolved())
                num_resolved = len(resolved)
                num_pending = len(var_entry.get_pending())
                if num_resolved == 1 and num_pending == 0:
                    # allele is resolved: yield and continue
                    a = var_entry.alleles[resolved[0]]
                    q = (
                        a.quality // a.length if a.length > 0 else 30
                    )  # Corner case empty ref allele
                    yield var_entry.variant_id, resolved[0], q
                elif num_resolved > 1 and num_pending == 0:
                    # multiple alleles possible: yield longest
                    lengths = [var_entry.alleles[r].length for r in resolved]
                    i = resolved[lengths.index(max(lengths))]
                    a = var_entry.alleles[i]
                    q = (
                        a.quality // a.length if a.length > 0 else 30
                    )  # Corner case empty ref allele
                    yield var_entry.variant_id, i, q
                elif num_pending > 0:
                    # allele is not resolved: re-queue
                    vqueue.appendleft(var_entry)
                    break
                # else: allele does match. discard and continue

        # After last cigar operation, yield ALL resolved variants, pop unresolved variants
        while vqueue:
            var_entry = vqueue.popleft()
            resolved = list(var_entry.get_resolved())
            num_resolved = len(resolved)
            num_pending = len(var_entry.get_pending())
            if num_resolved == 1 and num_pending == 0:
                # allele is resolved: yield and continue
                a = var_entry.alleles[resolved[0]]
                q = a.quality // a.length if a.length > 0 else 30  # Corner case empty ref allele
                yield var_entry.variant_id, resolved[0], q
            elif num_resolved > 1 and num_pending == 0:
                # multiple alleles possible: yield longest
                lengths = [var_entry.alleles[r].length for r in resolved]
                i = resolved[lengths.index(max(lengths))]
                a = var_entry.alleles[i]
                q = a.quality // a.length if a.length > 0 else 30  # Corner case empty ref allele
                yield var_entry.variant_id, i, q

    def detect_alleles_match(self, variant, entry, bam_read, ref_pos, query_pos, length):
        query_start = entry.query_start
        op_start = max(0, entry.query_start - query_pos)
        for i, a in enumerate(entry):
            # Skip already failed alleles
            if a.progress < 0:
                continue

            # Process remaining match-bases:
            ops_consumed = op_start
            allele_seq = variant.get_allele(i)
            query_pos = query_start + a.matched + a.inserted
            while a.matched < a.match_target and ops_consumed < length:
                qbase = bam_read.query_sequence[query_pos]
                vbase = allele_seq[a.matched + a.inserted]
                if qbase == vbase:
                    ops_consumed += 1
                    if bam_read.query_qualities:
                        a.quality += bam_read.query_qualities[query_pos]
                    else:
                        a.quality += 30  # TODO
                    a.matched += 1
                    a.progress += 1
                else:
                    break

            # If allele has non-matches left, but did not consume all match ops -> Fail allele
            if ops_consumed < length and a.progress < a.length:
                a.progress = -1

    def detect_alleles_insertion(self, variant, entry, bam_read, ref_pos, query_pos, length):
        query_start = entry.query_start
        for i, a in enumerate(entry):
            # Skip already failed alleles
            if a.progress < 0:
                continue

            # Process remaining insert ops:
            ops_consumed = 0
            allele_seq = variant.get_allele(i)
            while a.inserted < a.insert_target and ops_consumed < length:
                ops_consumed += 1
                qbase = bam_read.query_sequence[query_start + a.matched + a.inserted]
                vbase = allele_seq[a.matched + a.inserted]
                if qbase == vbase:
                    a.inserted += 1
                    a.progress += 1
                    a.quality += 30  # TODO
                else:
                    break

            # If allele has non-inserts left, but did not consume all insert ops -> Fail allele
            if ops_consumed < length and 0 < a.progress < a.length:
                a.progress = -1

    def detect_alleles_deletion(self, variant, entry, bam_read, ref_pos, query_pos, length):
        for i, a in enumerate(entry):
            # Skip already failed alleles
            if a.progress < 0:
                continue

            # Process remaining delete ops:
            ops_consumed = 0
            while a.deleted < a.delete_target and ops_consumed < length:
                ops_consumed += 1
                a.deleted += 1
                a.progress += 1
                a.quality += 30  # TODO

            # If allele has non-deletions left, but did not consume all delete ops -> Fail allele
            if ops_consumed < length and a.progress < a.length:
                a.progress = -1

    def detect_overlapping_variants(self, variants):
        """
        Checks for deletion variants overlapping other variants and for variants with duplicate
        positions. Returns a set of variant indices, which are conflict with another variant and
        should not be considered for allele detection.

        variants -- list of variants (VcfVariant objects)
        """
        j = 0
        conflicting = set()
        seen_pos = set()
        while j < len(variants):
            v = variants[j]
            if v.position in seen_pos:
                conflicting.add(j)
                j += 1
                continue
            else:
                seen_pos.add(v.position)
            ref = len(v.reference_allele)
            max_del = max(ref - len(alt) for alt in v.get_alt_allele_list())
            if max_del > 0:
                # at least one alt allele shorter than ref allele exists:
                deletion_end = v.position + ref
                if j + 1 < len(variants) and variants[j + 1].position < deletion_end:
                    # at least one follow-up variant overlaps the deletion
                    conflicting.add(j)
                    while j + 1 < len(variants) and variants[j + 1].position < deletion_end:
                        j += 1
                        conflicting.add(j)
            j += 1
        return conflicting

    def build_var_progress(self, variants, j):
        """
        Creates an object for tracking match progress of the j-th variant. Each object contains
        the variant id and lengths for every allele.
        """
        v = VariantProgress(j)
        ref_len = len(variants[j].reference_allele)
        v.add_allele(len(variants[j].reference_allele), 0, 0)
        for i, alt in enumerate(variants[j].get_alt_allele_list()):
            alt_len = len(alt)
            match_target = min(ref_len, alt_len)
            ins_target = max(0, len(alt) - ref_len)
            del_target = max(0, ref_len - len(alt))
            v.add_allele(match_target, ins_target, del_target)
        return v

    @staticmethod
    def split_cigar(cigar, i, consumed):
        """
        Split a CIGAR into two parts. i and consumed describe the split position.
        i is the element of the cigar list that should be split, and consumed says
        at how many operations to split within that element.

        The CIGAR is given as a list of (operation, length) pairs.

        i -- split at this index in cigar list
        consumed -- how many cigar ops at cigar[i] are to the *left* of the
            split position

        Return a tuple (left, right).

        Example:
        Assume the cigar is 3M 1D 6M 2I 4M.
        With i == 2 and consumed == 5, the cigar is split into
        3M 1D 5M and 1M 2I 4M.
        """
        middle_op, middle_length = cigar[i]
        assert consumed <= middle_length
        if consumed > 0:
            left = cigar[:i] + [(middle_op, consumed)]
        else:
            left = cigar[:i]
        if consumed < middle_length:
            right = [(middle_op, middle_length - consumed)] + cigar[i + 1 :]
        else:
            right = cigar[i + 1 :]
        return left, right

    @staticmethod
    def cigar_prefix_length(cigar, reference_bases):
        """
        Given a prefix of length reference_bases relative to the reference, how
        long is the prefix of the read? In other words: If reference_bases on
        the reference are consumed, how many bases on the query does that
        correspond to?

        If the position is within or at the end of an insertion (which do not
        consume bases on the reference), then the number of bases up to the
        beginning of the insertion is reported.

        Return a pair (reference_bases, query_bases) where the value for
        reference_bases may be smaller than the requested one if the CIGAR does
        not cover enough reference bases.

        Reference skips (N operators) are treated as the end of the read. That
        is, no positions beyond a reference skip are reported.
        """
        ref_pos = 0
        query_pos = 0
        for op, length in cigar:
            if op in (0, 7, 8):  # M, X, =
                ref_pos += length
                query_pos += length
                if ref_pos >= reference_bases:
                    return (reference_bases, query_pos + reference_bases - ref_pos)
            elif op == 2:  # D
                ref_pos += length
                if ref_pos >= reference_bases:
                    return (reference_bases, query_pos)
            elif op == 1:  # I
                query_pos += length
            elif op == 4 or op == 5:  # soft or hard clipping
                pass
            elif op == 3:  # N
                # Always stop at reference skips
                return (reference_bases, query_pos)
            else:
                assert False, "unknown CIGAR operator"
        assert ref_pos < reference_bases
        return (ref_pos, query_pos)

    @staticmethod
    def realign(
        variant,
        bam_read,
        cigartuples,
        i,
        consumed,
        query_pos,
        reference,
        overhang,
        use_affine,
        gap_start,
        gap_extend,
        default_mismatch,
    ):
        """
        Realign a read to the two alleles of a single variant.
        i and consumed describe where to split the cigar into a part before the
        variant position and into a part starting at the variant position, see split_cigar().

        variant -- VcfVariant
        bam_read -- the AlignedSegment
        cigartuples -- the AlignedSegment.cigartuples property (accessing it is expensive, so re-use it)
        i, consumed -- see split_cigar method
        query_pos -- index of the query base that is at the variant position
        reference -- the reference as a str-like object (full chromosome)
        overhang -- extend alignment by this many bases to left and right
        use_affine -- if true, use affine gap costs for realignment
        gap_start, gap_extend -- if affine_gap=true, use these parameters for affine gap cost alignment
        default_mismatch -- if affine_gap=true, use this as mismatch cost in case no base qualities are in bam
        """
        # Do not process symbolic alleles like <DEL>, <DUP>, etc.
        if variant.alternative_allele.startswith("<"):
            return None, None

        left_cigar, right_cigar = ReadSetReader.split_cigar(cigartuples, i, consumed)

        left_ref_bases, left_query_bases = ReadSetReader.cigar_prefix_length(
            left_cigar[::-1], overhang
        )
        right_ref_bases, right_query_bases = ReadSetReader.cigar_prefix_length(
            right_cigar, len(variant.reference_allele) + overhang
        )

        assert variant.position - left_ref_bases >= 0
        assert variant.position + right_ref_bases <= len(reference)

        query = bam_read.query_sequence[
            query_pos - left_query_bases : query_pos + right_query_bases
        ]
        ref = reference[variant.position - left_ref_bases : variant.position + right_ref_bases]
        alt = (
            reference[variant.position - left_ref_bases : variant.position]
            + variant.alternative_allele
            + reference[
                variant.position
                + len(variant.reference_allele) : variant.position
                + right_ref_bases
            ]
        )

        if use_affine:
            assert gap_start is not None
            assert gap_extend is not None
            assert default_mismatch is not None

            # get base qualities if present (to be used as mismatch costs)
            base_qualities = [default_mismatch] * len(query)
            # if bam_read.query_qualities != None:
            #    base_qualities = bam_read.query_qualities[query_pos-left_query_bases:query_pos+right_query_bases]

            # compute edit dist. with affine gap costs using base qual. as mismatch cost
            distance_ref = edit_distance_affine_gap(
                query, ref, base_qualities, gap_start, gap_extend
            )
            distance_alt = edit_distance_affine_gap(
                query, alt, base_qualities, gap_start, gap_extend
            )
            base_qual_score = abs(distance_ref - distance_alt)
        else:
            base_qual_score = 30
            distance_ref = edit_distance(query, ref)
            distance_alt = edit_distance(query, alt)

        if distance_ref < distance_alt:
            return 0, base_qual_score  # detected REF
        elif distance_ref > distance_alt:
            return 1, base_qual_score  # detected ALT
        else:
            return None, None  # cannot decide

    @staticmethod
    def detect_alleles_by_alignment(
        variants,
        j,
        bam_read,
        reference,
        overhang=10,
        use_affine=False,
        gap_start=None,
        gap_extend=None,
        default_mismatch=None,
    ):
        """
        Detect which alleles the given bam_read covers. Detect the correct
        alleles of the variants that are covered by the given bam_read.

        Yield tuples (position, allele, quality).

        variants -- list of variants (VcfVariant objects)
        j -- index of the first variant (in the variants list) to check
        """
        # Accessing bam_read.cigartuples is expensive, do it only once
        cigartuples = bam_read.cigartuples

        # For the same reason, the following check is here instad of
        # in the _usable_alignments method
        if not cigartuples:
            return

        for index, i, consumed, query_pos in _iterate_cigar(variants, j, bam_read, cigartuples):
            allele, quality = ReadSetReader.realign(
                variants[index],
                bam_read,
                cigartuples,
                i,
                consumed,
                query_pos,
                reference,
                overhang,
                use_affine,
                gap_start,
                gap_extend,
                default_mismatch,
            )
            if allele in (0, 1):
                yield (index, allele, quality)  # TODO quality???

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        self._reader.close()


def merge_two_reads(read1: Read, read2: Read) -> Read:
    """
    Merge two reads *that belong to the same haplotype* (such as the two
    ends of a paired-end read) into a single Read. Overlaps are allowed.
    """
    assert read1.is_sorted()
    assert read2.is_sorted()
    if read2:
        result = Read(
            read1.name,
            read1.mapqs[0],
            read1.source_id,
            read1.sample_id,
            read1.reference_start,
            read1.BX_tag,
        )
        result.add_mapq(read2.mapqs[0])
    else:
        return read1

    i1 = 0
    i2 = 0

    def add1():
        result.add_variant(read1[i1].position, read1[i1].allele, read1[i1].quality)

    def add2():
        result.add_variant(read2[i2].position, read2[i2].allele, read2[i2].quality)

    while i1 < len(read1) or i2 < len(read2):
        if i1 == len(read1):
            add2()
            i2 += 1
            continue
        if i2 == len(read2):
            add1()
            i1 += 1
            continue
        variant1 = read1[i1]
        variant2 = read2[i2]
        if variant2.position < variant1.position:
            add2()
            i2 += 1
        elif variant2.position > variant1.position:
            add1()
            i1 += 1
        else:
            # Variant on self-overlapping read pair
            assert read1[i1].position == read2[i2].position
            # If both alleles agree, merge into single variant and add up qualities
            if read1[i1].allele == read2[i2].allele:
                quality = read1[i1].quality + read2[i2].quality
                result.add_variant(read1[i1].position, read1[i1].allele, quality)
            else:
                # Otherwise, take variant with highest base quality and discard the other.
                if read1[i1].quality >= read2[i2].quality:
                    add1()
                else:
                    add2()
            i1 += 1
            i2 += 1
    return result


def merge_reads(*reads: Read) -> Read:
    """
    Merge multiple reads that belong to the same haplotype into a single Read.

    If the iterable is empty, a ValueError is raised.

    This 'naive' version just calls merge_two_reads repeatedly on all the reads.

    # TODO
    # The actual challenge is dealing with conflicts in variants covered by
    # more than one read. A solution would be to not merge if there are any
    # (or too many) conflicts and let the main algorithm deal with it.
    """
    it = iter(reads)
    try:
        read = next(it)
    except StopIteration:
        raise ValueError("no reads to merge")
    assert read.is_sorted()
    for partner in it:
        read = merge_two_reads(read, partner)
    return read
