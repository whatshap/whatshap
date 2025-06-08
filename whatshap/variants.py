"""
Detect variants in reads.
"""

import logging
import csv
from collections import defaultdict, Counter
from typing import Iterable, Iterator, List, Optional
from dataclasses import dataclass

import pysam
from pysam import AlignedSegment

from .vcf import VcfVariant
from .bam import SampleBamReader, MultiBamReader, BamReader, AlignmentWithSourceID
from .core import Genotype, Read, ReadSet, NumericSampleIds
from .align import edit_distance, edit_distance_affine_gap, kmer_align, enumerate_all_kmers
from ._variants import _iterate_cigar, _detect_alleles

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


@dataclass
class AlignedRead:
    read: Read
    is_supplementary: bool
    is_reverse: bool
    reference_start: int
    reference_end: int

    def distance(self, other: "AlignedRead") -> int:
        return max(
            0,
            self.reference_start - other.reference_end,
            other.reference_start - self.reference_end,
        )


def is_alignment_primary(alignment: pysam.AlignedSegment) -> bool:
    non_primary: bool = False
    non_primary |= alignment.is_supplementary
    non_primary |= alignment.is_secondary
    non_primary |= alignment.is_unmapped
    return not non_primary


def is_alignmentwsid_primary(alignment: AlignmentWithSourceID) -> bool:
    return is_alignment_primary(alignment=alignment.bam_alignment)


# some value that would very unlikely appear as a suffix of a genuine read id
PRIMARY_DEFAULT_SUB_ALIGNMENT_ID = "____1"


def get_sub_alignment_id(
    alignment: pysam.AlignedSegment,
    is_primary: bool,
    primary_default: Optional[str] = PRIMARY_DEFAULT_SUB_ALIGNMENT_ID,
) -> str:
    """
    For every alignment segment we want to have a unique alignment segment id, while allowing for a "guessable"
     id for a primary alignment. The primary point here is to the have ability to differentiate between supplementary
     alignments for a given read. CIGAR string provides, almost always, a unique differentiation between supplementary
     alignments, except for an identical subalignment where the start/end of the reads aligns identically to the same
     region of the reference. For such a case the flag is also considered for the id generation.
    We are using the string representation vs int ones for cpp/python interop, as the python hash int values are too
     large to handle in standard int type in cpp.
    """
    if not is_primary or primary_default is None:
        return str(hash((alignment.cigarstring, alignment.flag)))
    return primary_default


def get_sub_alignmentw_id_wsid(
    alignment: AlignmentWithSourceID,
    is_primary: bool,
    primary_default: Optional[str] = PRIMARY_DEFAULT_SUB_ALIGNMENT_ID,
) -> str:
    return get_sub_alignment_id(
        alignment=alignment.bam_alignment, is_primary=is_primary, primary_default=primary_default
    )


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
        *,
        mapq_threshold: int = 20,
        overhang: int = 10,
        affine: int = False,
        gap_start: int = 10,
        gap_extend: int = 7,
        default_mismatch: int = 15,
        duplicates: bool = False,
        use_kmerald: bool = False,
        kmeralign_costs_path: Optional[str] = None,
        kmer_size: int = 7,
        kmerald_gappenalty: float = 40,
        kmerald_window: int = 25,
        use_supplementary: bool = False,
        supplementary_distance_threshold: int = 100_000,
        allow_supplementary_only_read_groups: bool = False,
    ):
        """
        paths -- list of BAM paths
        reference -- path to reference FASTA (can be None)
        numeric_sample_ids -- ??
        mapq_threshold -- minimum mapping quality
        overhang -- extend alignment by this many bases to left and right
        affine -- use affine gap costs
        gap_start, gap_extend, default_mismatch -- parameters for affine gap cost alignment
        duplicates -- read alignments marked as duplicate
        allow_supplementary_only_read_groups - if we are allowing for independent supplementary alignments consideration
        """
        self._mapq_threshold = mapq_threshold
        self._numeric_sample_ids = numeric_sample_ids
        self._use_affine = affine
        self._gap_start = gap_start
        self._gap_extend = gap_extend
        self._default_mismatch = default_mismatch
        self._overhang = overhang
        self._duplicates = duplicates
        self._use_kmerald = use_kmerald
        self._kmeralign_costs_path = kmeralign_costs_path
        self._kmer_size = kmer_size
        self._kmerald_gappenalty = kmerald_gappenalty
        self._kmerald_window = kmerald_window
        self._paths = paths
        self._reader: BamReader
        self._use_supplementary = use_supplementary
        self._supplementary_distance_threshold = supplementary_distance_threshold
        self._allow_supplementary_only_read_groups = allow_supplementary_only_read_groups
        if len(paths) == 1:
            self._reader = SampleBamReader(paths[0], reference=reference)
        else:
            self._reader = MultiBamReader(paths, reference=reference)

    @property
    def n_paths(self) -> int:
        return len(self._paths)

    def read(
        self,
        chromosome,
        variants,
        sample,
        reference,
        regions=None,
        restricted_genotypes: Optional[List[Genotype]] = None,
    ) -> ReadSet:
        """
        Detect alleles and return a ReadSet object containing reads representing
        the given variants.

        If a reference is provided (reference is not None), alleles are
        detected by re-aligning sections of the query to the REF and ALT
        sequence extended a few bases to the left and right.

        If reference is None, alleles are detected by inspecting the
        existing alignment (via the CIGAR).

        If the restricted_genotypes is not None, then only alleles from restricted_genotypes[i] will be considered for variant i. The length of restricted_genotypes should match the length of variants.

        chromosome -- name of chromosome to work on
        variants -- list of vcf.VcfVariant objects
        sample -- name of sample to work on. If None, read group information is
            ignored and all reads in the file are used.
        reference -- reference sequence of the given chromosome (or None)
        regions -- list of start,end tuples (end can be None)
        restricted_genotypes -- list of restricted genotypes (or None if there is no reliable auxiliary information).
        """
        # Since variants are identified by position, positions must be unique.
        if __debug__ and variants:
            varposc = Counter(variant.position for variant in variants)
            pos, count = varposc.most_common()[0]
            assert count == 1, f"Position {pos} occurs more than once in variant list."

        assert restricted_genotypes is None or len(restricted_genotypes) == len(variants)
        alignments = self._usable_alignments(chromosome, sample, regions)
        reads = self._alignments_to_reads(
            alignments, variants, sample, reference, restricted_genotypes
        )
        grouped_reads = self._group_reads(
            reads,
            self._supplementary_distance_threshold,
            allow_supplementary_only_groups=self._allow_supplementary_only_read_groups,
        )
        readset = self._make_readset_from_grouped_reads(grouped_reads)
        return readset

    @staticmethod
    def _make_readset_from_grouped_reads(groups: Iterable[List[Read]]) -> ReadSet:
        read_set = ReadSet()
        for group in groups:
            read_set.add(merge_reads(*group))
        return read_set

    @staticmethod
    def create_read_from_group(
        group: List[AlignedRead],
        distance_threshold: int,
        allow_supplementary_only_groups: bool = False,
    ) -> Optional[Read]:
        """
        Merge multiple AlignedReads into a single Read.

        Pick supplementary reads that have the same orientation as the primary and with
        the distance at most distance_threshold from primary, find the set of variants that fully agree
        and return these variants as a read.

        If the list does not contain primary reads, then return None, unless supplementary only groups are allowed.
        If the list contains more than two primary alignments, return None and report a warning.
        """
        if len(group) > 1:
            logger.debug(f"Group of read {group[0].read.name!r} has {len(group)} items.")
        primary: Optional[AlignedRead] = None
        n_primary = 0
        for read in group:
            if not read.is_supplementary:
                n_primary += 1
                primary = read
        if primary is None and not allow_supplementary_only_groups:
            return None
        if n_primary > 2:
            logger.warning(
                f"Read name {group[0].read.name!r} has more than two primary alignments."
            )
            return None
        if primary is None:
            # we can only be here is we allow for supplementary only groups and there was no primary,
            #  in that scenario we would expect a single-read group input
            #  so we can assign any, i.e., the supplementary read in the group as a "primary"
            primary = read
        reference_start = primary.reference_start
        variants = dict()
        skip = set()
        for read in group:
            if read.is_supplementary:
                if read.is_reverse != primary.is_reverse:
                    continue
                if primary.distance(read) > distance_threshold:
                    continue
            reference_start = min(reference_start, read.reference_start)
            for variant in read.read:
                if variant.position in variants:
                    if variants[variant.position].allele != variant.allele:
                        skip.add(variant.position)
                else:
                    variants[variant.position] = variant
        union_read = Read(
            # because of the way the ReadSet datastructure is utilized, we need to ensure that we can
            #   store not only primary but supplementary read alignments in the same ReadSet datastructure
            primary.read.name
            + (primary.read.sub_alignment_id if allow_supplementary_only_groups else ""),
            primary.read.mapqs[0],
            primary.read.source_id,
            primary.read.sample_id,
            reference_start,
            primary.read.BX_tag,
            primary.read.HP_tag,
            primary.read.PS_tag,
            chromosome=primary.read.chromosome,
            sub_alignment_id=primary.read.sub_alignment_id,
            is_supplementary=primary.read.is_supplementary,
            is_reverse=primary.is_reverse,
            reference_end=primary.reference_end,
        )
        for position, variant in variants.items():
            if position not in skip:
                union_read.add_variant(variant.position, variant.allele, variant.quality)
        union_read.sort()
        if len(union_read) != len(primary.read):
            logger.debug(
                f"Converted read {primary.read.name} with {len(primary.read)} variants"
                f" to read with {len(union_read)} variants."
            )
        return union_read

    @staticmethod
    def _group_reads(
        reads: Iterable[AlignedRead],
        distance_threshold: int,
        allow_supplementary_only_groups: bool = False,
    ) -> Iterator[List[Read]]:
        """
        Group reads into paired-end read pairs. Uses name, source_id and sample_id
        as grouping key.

        TODO
        Grouping by name should be sufficient since the SAM spec states:
        "Reads/segments having identical QNAME are regarded to come from the same template."
        """
        groups = defaultdict(list)
        for read in reads:
            groups[
                (
                    read.read.source_id,
                    read.read.name,
                    # by using this sub-alignment-id trick
                    # when allowing for supplementary only groups
                    # we are guaranteed to have single-read only groups
                    read.read.sub_alignment_id if allow_supplementary_only_groups else None,
                    read.read.sample_id,
                )
            ].append(read)
        n_skipped = 0
        n_non_singleton = 0
        for group in groups.values():
            if len(group) > 1:
                n_non_singleton += 1
            read = ReadSetReader.create_read_from_group(
                group,
                distance_threshold,
                allow_supplementary_only_groups=allow_supplementary_only_groups,
            )
            if read is None:
                n_skipped += 1
            else:
                yield [read]

        logger.info(f"Number of non-singleton groups: {n_non_singleton}")
        logger.info(f"Skipped {n_skipped} groups")

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
                    (not self._use_supplementary and alignment.bam_alignment.is_supplementary)
                    or alignment.bam_alignment.mapping_quality < self._mapq_threshold
                    or alignment.bam_alignment.is_secondary
                    or alignment.bam_alignment.is_unmapped
                    or (not self._duplicates and alignment.bam_alignment.is_duplicate)
                ):
                    continue
                yield alignment

    def has_reference(self, chromosome):
        return self._reader.has_reference(chromosome)

    def _alignments_to_reads(
        self,
        alignments,
        variants,
        sample,
        reference,
        restricted_genotypes: Optional[List[Genotype]],
    ):
        """
        Convert BAM alignments to Read objects.

        If reference is not None, alleles are detected through re-alignment.

        Yield Read objects.
        """
        # FIXME hard-coded zero
        numeric_sample_id = 0 if sample is None else self._numeric_sample_ids[sample]
        number_of_alignments = 0
        number_of_supplementary_alignments = 0
        if reference is not None:
            # Copy the pyfaidx.FastaRecord into a str for faster access
            reference = reference[:]
            normalized_variants = variants
        else:
            normalized_variants = [variant.normalized() for variant in variants]

        # Create allele progress trackers once, instead of doing it for every read again
        if reference is None:
            # Discard overlapping and duplicate-positioned variants for more efficient iteration
            valid_variant_ids = self.detect_non_overlapping_variants(normalized_variants)
            valid_positions = [normalized_variants[j].position for j in valid_variant_ids]
            var_progress = [
                self.build_var_progress(normalized_variants, j) for j in valid_variant_ids
            ]
            var_progress.sort(key=lambda x: x.variant_id)

        i = 0  # index into variants (reference) or variant progresses (no reference)

        if self._use_kmerald:
            calculated_costs = {}
            splitted_strings = {}
            kmerald_costs = {}
            with open(self._kmeralign_costs_path) as costs_file:
                reader = csv.reader(costs_file, delimiter="\t")
                for line in reader:
                    kmerald_costs[(int(line[0]), int(line[1]))] = line[2]
        else:
            kmerald_costs = None
            calculated_costs = None
            splitted_strings = None

        def get_tag_or_default(aln, tag, default):
            if aln.bam_alignment.has_tag(tag):
                return aln.bam_alignment.get_tag(tag)
            else:
                return default

        for alignment in alignments:
            barcode = get_tag_or_default(alignment, "BX", "")
            hp = get_tag_or_default(alignment, "HP", -1)
            ps = get_tag_or_default(alignment, "PS", -1)
            # Some 10X bams appear to use a Z type for the PS tag even when they are integers
            # so here we try to cast to an integer, but if that fails we raise an error
            try:
                ps = int(ps)
            except ValueError:
                raise ValueError(
                    f"Invalid PS tag value ({ps}) in read {alignment.bam_alignment.query_name}. PS must be an integer."
                )
            primary: bool = is_alignmentwsid_primary(alignment)
            sub_alignment_id: str = get_sub_alignment_id(
                alignment.bam_alignment, is_primary=primary
            )
            read = Read(
                alignment.bam_alignment.query_name,
                alignment.bam_alignment.mapq,
                alignment.source_id,
                numeric_sample_id,
                alignment.bam_alignment.reference_start,
                barcode,
                hp,
                ps,
                chromosome=alignment.bam_alignment.reference_name,
                sub_alignment_id=sub_alignment_id,
                is_supplementary=alignment.bam_alignment.is_supplementary,
                is_reverse=alignment.bam_alignment.is_reverse,
                reference_end=alignment.bam_alignment.reference_end,
            )

            if reference is None:
                # Skip variant progress objects that are to the left of this read
                while (
                    i < len(valid_positions)
                    and valid_positions[i] < alignment.bam_alignment.reference_start
                ):
                    i += 1
                detected = _detect_alleles(
                    normalized_variants, var_progress, i, alignment.bam_alignment
                )
            else:
                # Skip variants that are to the left of this read
                while (
                    i < len(normalized_variants)
                    and normalized_variants[i].position < alignment.bam_alignment.reference_start
                ):
                    i += 1
                detected = self.detect_alleles_by_alignment(
                    variants,
                    restricted_genotypes,
                    i,
                    alignment.bam_alignment,
                    reference,
                    self._overhang,
                    self._use_affine,
                    self._gap_start,
                    self._gap_extend,
                    self._default_mismatch,
                    self._use_kmerald,
                    kmerald_costs,
                    self._kmer_size,
                    self._kmerald_gappenalty,
                    self._kmerald_window,
                    calculated_costs,
                    splitted_strings,
                )

            for j, allele, quality in detected:
                read.add_variant(variants[j].position, allele, quality)
            if read:  # At least one variant covered and detected
                number_of_alignments += 1
                number_of_supplementary_alignments += alignment.bam_alignment.is_supplementary
                yield AlignedRead(
                    read,
                    alignment.bam_alignment.is_supplementary,
                    alignment.bam_alignment.is_reverse,
                    alignment.bam_alignment.reference_start,
                    alignment.bam_alignment.reference_end,
                )

        logger.info(f"Number of supplementary alignments: {number_of_supplementary_alignments}")

    def detect_non_overlapping_variants(self, variants: List[VcfVariant]):
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
        return [j for j in range(len(variants)) if j not in conflicting]

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
    def split_cigar_left(cigar, i, consumed):
        """
        Split a CIGAR into two parts. i and consumed describe the split position.
        i is the element of the cigar list that should be split, and consumed says
        at how many operations to split within that element.

        The CIGAR is given as a list of (operation, length) pairs.

        i -- split at this index in cigar list
        consumed -- how many cigar ops at cigar[i] are to the *left* of the
            split position

        Example:
        Assume the cigar is 3M 1D 6M 2I 4M.
        With i == 2 and consumed == 5, the cigar is split into
        3M 1D 5M and 1M 2I 4M.

        This function returns the left part of the split as an iterator that runs
        from the split point to the beginning of the full cigar.
        See split_cigar_right to obtain the right part.
        """
        middle_op, middle_length = cigar[i]
        assert consumed <= middle_length
        if consumed > 0:
            yield middle_op, consumed
        for j in range(i - 1, -1, -1):
            yield cigar[j]

    @staticmethod
    def split_cigar_right(cigar, i, consumed):
        """
        Counterpart of split_cigar_left that returns an iterator for the right
        part of the split, running from the split point to the end of the input
        cigar.
        """
        middle_op, middle_length = cigar[i]
        if consumed < middle_length:
            yield middle_op, middle_length - consumed
        for j in range(i + 1, len(cigar)):
            yield cigar[j]

    @staticmethod
    def cigar_prefix_length(cigar, reference_bases: int):
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
        variant: VcfVariant,
        restricted_variants: Optional[Genotype],
        bam_read: AlignedSegment,
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
        use_kmerald,
        kmerald_costs,
        kmer_size,
        kmerald_gappenalty,
        kmerald_window,
        calculated_costs,
        splitted_strings,
    ):
        """
        Realign a read to the two alleles (or to the alleles from genotype if it is not None) of a single variant.
        i and consumed describe where to split the cigar into a part before the
        variant position and into a part starting at the variant position, see split_cigar().

        variant -- VcfVariant
        restricted_variants -- genotype with restricted variants (or None if there is no such information)
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
        if any(alt.startswith("<") for alt in variant.get_alt_allele_list()):
            return None, None

        left_cigar_iterator = ReadSetReader.split_cigar_left(cigartuples, i, consumed)
        right_cigar_iterator = ReadSetReader.split_cigar_right(cigartuples, i, consumed)

        if use_kmerald:
            left_ref_bases, left_query_bases = ReadSetReader.cigar_prefix_length(
                left_cigar_iterator, int(kmerald_window)
            )
            right_ref_bases, right_query_bases = ReadSetReader.cigar_prefix_length(
                right_cigar_iterator, len(variant.reference_allele) + int(kmerald_window)
            )
            assert variant.position - left_ref_bases >= 0
            assert variant.position + right_ref_bases <= len(reference)
            query_temp = bam_read.query_sequence[
                query_pos - left_query_bases : query_pos + right_query_bases
            ]
            if query_temp in splitted_strings:
                query = splitted_strings[query_temp]
            else:
                query = enumerate_all_kmers(str(query_temp).encode("UTF-8"), int(kmer_size))
                splitted_strings[query_temp] = query

            ref_temp = reference[
                variant.position - left_ref_bases : variant.position + right_ref_bases
            ]
            if ref_temp in splitted_strings:
                ref = splitted_strings[ref_temp]
            else:
                ref = enumerate_all_kmers(str(ref_temp).encode("UTF-8"), int(kmer_size))
                splitted_strings[ref_temp] = ref

            alt_temp = (
                reference[variant.position - left_ref_bases : variant.position]
                + variant.alternative_allele
                + reference[
                    variant.position
                    + len(variant.reference_allele) : variant.position
                    + right_ref_bases
                ]
            )

            if alt_temp in splitted_strings:
                alt = splitted_strings[alt_temp]
            else:
                alt = enumerate_all_kmers(str(alt_temp).encode("UTF-8"), int(kmer_size))
                splitted_strings[alt_temp] = alt

            base_qual_score = 30
            distance_ref = 0
            distance_alt = 0
            if (ref_temp, query_temp) in calculated_costs:
                distance_ref = calculated_costs[(ref_temp, query_temp)]
            else:
                distance_ref = kmer_align(ref, query, kmerald_costs, kmerald_gappenalty)
                calculated_costs[(ref_temp, query_temp)] = distance_ref

            if (alt_temp, query_temp) in calculated_costs:
                distance_alt = calculated_costs[(alt_temp, query_temp)]
            else:
                distance_alt = kmer_align(alt, query, kmerald_costs, kmerald_gappenalty)
                calculated_costs[(alt_temp, query_temp)] = distance_alt

            if distance_ref < distance_alt:
                return 0, base_qual_score  # detected REF
            elif distance_ref > distance_alt:
                return 1, base_qual_score  # detected ALT
            else:
                return None, None  # cannot decide
        else:
            left_ref_bases, left_query_bases = ReadSetReader.cigar_prefix_length(
                left_cigar_iterator, overhang
            )
            right_ref_bases, right_query_bases = ReadSetReader.cigar_prefix_length(
                right_cigar_iterator, len(variant.reference_allele) + overhang
            )

            assert variant.position - left_ref_bases >= 0
            assert variant.position + right_ref_bases <= len(reference)

            query = bam_read.query_sequence[
                query_pos - left_query_bases : query_pos + right_query_bases
            ]
            pos = variant.position
            left_pad = reference[pos - left_ref_bases : pos]
            right_pad = reference[pos + len(variant.reference_allele) : pos + right_ref_bases]
            padded_alleles = [reference[pos - left_ref_bases : pos + right_ref_bases]]
            for alt in variant.get_alt_allele_list():
                padded_alleles.append(left_pad + alt + right_pad)

        if use_affine:
            assert gap_start is not None
            assert gap_extend is not None
            assert default_mismatch is not None

            # get base qualities if present (to be used as mismatch costs)
            base_qualities = [default_mismatch] * len(query)
            # if bam_read.query_qualities != None:
            #    base_qualities = bam_read.query_qualities[query_pos-left_query_bases:query_pos+right_query_bases]

            # compute edit dist. with affine gap costs using base qual. as mismatch cost
            distances = [
                (i, edit_distance_affine_gap(query, allele, base_qualities, gap_start, gap_extend))
                for i, allele in enumerate(padded_alleles)
                if restricted_variants is None or i in restricted_variants.as_vector()
            ]
            distances.sort(key=lambda x: x[1])
            base_qual_score = (
                distances[0][1] - distances[1][1] if len(distances) > 1 else distances[0][1]
            )
        else:
            distances = [
                (i, edit_distance(query, allele))
                for i, allele in enumerate(padded_alleles)
                if restricted_variants is None or i in restricted_variants.as_vector()
            ]
            distances.sort(key=lambda x: x[1])
            base_qual_score = 30

        if len(distances) == 1 or distances[0][1] < distances[1][1]:
            return distances[0][0], base_qual_score  # detected REF
        else:
            return None, None  # cannot decide

    @staticmethod
    def detect_alleles_by_alignment(
        variants: List[VcfVariant],
        restricted_genotypes: Optional[List[Genotype]],
        j,
        bam_read: AlignedSegment,
        reference,
        overhang=10,
        use_affine=False,
        gap_start=None,
        gap_extend=None,
        default_mismatch=None,
        use_kmerald=False,
        kmerald_costs=None,
        kmer_size=7,
        kmerald_gappenalty=40,
        kmerald_window=25,
        calculated_costs=None,
        splitted_strings=None,
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

        # For the same reason, the following check is here instead of
        # in the _usable_alignments method
        if not cigartuples:
            return

        for index, i, consumed, query_pos in _iterate_cigar(variants, j, bam_read, cigartuples):
            allele, quality = ReadSetReader.realign(
                variants[index],
                restricted_genotypes[index] if restricted_genotypes else None,
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
                use_kmerald,
                kmerald_costs,
                kmer_size,
                kmerald_gappenalty,
                kmerald_window,
                calculated_costs,
                splitted_strings,
            )
            num_alts = len(variants[index].get_alt_allele_list())
            if allele in range(num_alts + 1):
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
            read1.HP_tag,
            read1.PS_tag,
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
