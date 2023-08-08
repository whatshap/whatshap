# cython: language_level=3

import logging
from collections import deque


logger = logging.getLogger(__name__)


def _iterate_cigar(variants, int j, bam_read, cigartuples):
    """
    Iterate over the CIGAR of the given bam_read and variants[j:] in lockstep.

    Yield tuples (index, i, consumed, query_pos) where index is into the variants list

    i and consumed describe the split position in the cigar

    bam_read -- a pysam.AlignedSegment
    variants -- list of variants (VcfVariant objects)
    j -- index of the first variant (in the variants list) to check
    """
    cdef:
        int ref_pos = bam_read.reference_start  # position relative to reference
        int query_pos = 0  # position relative to read
        int cigar_op
        int length
        int i
        int n = len(variants)
        int v_position

    # Skip variants that are located to the left of the read
    while j < n and variants[j].position < ref_pos:
        j += 1

    # Iterate over the CIGAR sequence (defining the alignment) and variant list in lockstep
    for i, (cigar_op, length) in enumerate(cigartuples):
        # The mapping of CIGAR operators to numbers is:
        # MIDNSHPX= => 012345678
        if j < n:
            v_position = variants[j].position
        if cigar_op in (0, 7, 8):  # M, X, = operators (match)
            # Iterate over all variants that are in this matching region
            while j < n and v_position < ref_pos + length:
                assert v_position >= ref_pos
                yield (j, i, v_position - ref_pos, query_pos + v_position - ref_pos)
                j += 1
                if j < n:
                    v_position = variants[j].position
            query_pos += length
            ref_pos += length
        elif cigar_op == 1:  # I operator (insertion)
            # TODO it should work to *not* handle the variant here, but at the next M or D region
            if j < n and v_position == ref_pos:
                yield (j, i, 0, query_pos)
                j += 1
                if j < n:
                    v_position = variants[j].position
            query_pos += length
        elif cigar_op == 2:  # D operator (deletion)
            # Iterate over all variants that are in this deleted region
            while j < n and v_position < ref_pos + length:
                assert v_position >= ref_pos
                yield (j, i, v_position - ref_pos, query_pos)
                j += 1
                if j < n:
                    v_position = variants[j].position
            ref_pos += length
        elif cigar_op == 3:  # N operator (reference skip)
            # Iterate over all variants that are in this skipped region
            while j < n and v_position < ref_pos + length:
                assert v_position >= ref_pos
                j += 1
                if j < n:
                    v_position = variants[j].position
            ref_pos += length
        elif cigar_op == 4:  # S operator (soft clipping)
            query_pos += length
        elif cigar_op == 5 or cigar_op == 6:  # H or P (hard clipping or padding)
            pass
        else:
            raise ValueError("Unsupported CIGAR operation: {}".format(cigar_op))


def _detect_alleles(variants, var_progress, first, bam_read):
    """
    Detect the correct alleles of the variants that are covered by the
    given bam_read.

    Yield tuples (index, allele, quality), where index is into the variants list.

    variants -- list of variants (VcfVariant objects)
    var_progress -- list of VariantProgress objects, also representing non-conflict variant ids
    first -- index of the first variant (in the var_progress list) to check
    bam_read -- alignment object to query sequence and qualities
    """
    cdef:
        int ref_pos = bam_read.reference_start  # position relative to reference
        int query_pos = 0                       # position relative to read
        int ref_end                             # end of ref span for cigar op
        int query_end                           # end of query span for cigar op
        int ref_len
        int j = first                           # index into var_progress
        int var_id                              # index into variants
        int var_pos                             # position of var_id-th variant
        int n = len(var_progress)
        int cigar_op                            # copy python vars here ...
        int length                              # ... for runtime optimization

    # Skip variants that come before this region
    while j < n:
        var_id = var_progress[j].variant_id
        var_pos = variants[var_id].position
        if var_pos >= ref_pos:
            break
        j += 1

    vqueue = deque()  # buffer for pending variants to keep them in positional order

    for py_cigar_op, py_length in bam_read.cigartuples:
        cigar_op, length = py_cigar_op, py_length  # much faster when using typed cython variables

        # Skip variants that come before this region
        while j < n:
            var_id = var_progress[j].variant_id
            var_pos = variants[var_id].position
            if var_pos >= ref_pos:
                break
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
        while j < n:
            var_id = var_progress[j].variant_id
            var_pos = variants[var_id].position
            # Stop when exceeding end of cigar op span
            if var_pos >= ref_end:
                break

            ref_len = len(variants[var_id].reference_allele)
            # Special case: If a non-insertion variant is seen by I-Op, continue with next Op
            # It is an insertion in front of a non-insertion variant, that must be ignored
            # We cannot skip I-Op in general, because this might overlook insertion variants
            if cigar_op == 1 and ref_len > 0:
                break
            # Special case: If a D-Op sees an insertion variant, skip this variant to be conform
            # with old implementation. Actually, it would be correct to assume ref allele here,
            # if the preivous base matched. This seems to be an artifact of normalized variants.
            if cigar_op == 2 and ref_len == 0:
                j += 1
                continue

            # Entry = [var_id, query_start_pos, for each allele: AlleleProgress object]
            query_start = query_pos + var_pos - ref_pos if cigar_op != 2 else query_pos
            var_progress[j].reset(query_start)
            vqueue.append(var_progress[j])
            j += 1

        # Handle detection and positional progress depending on cigar op
        ref_end = ref_pos
        query_end = query_pos
        if cigar_op in (0, 7, 8):  # M, X, = operators (match)
            handler = _detect_alleles_match
            ref_end += length
            query_end += length
        elif cigar_op == 1:  # I operator (insertion)
            handler = _detect_alleles_insertion
            query_end += length
        elif cigar_op == 2:  # D operator (deletion)
            handler = _detect_alleles_deletion
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
            if num_resolved >= 1 and num_pending == 0:
                # multiple alleles possible: yield longest
                lengths = [var_entry.alleles[r].length for r in resolved]
                i = resolved[lengths.index(max(lengths))]
                a = var_entry.alleles[i]
                q = a.quality // a.length if a.length > 0 else 30  # Corner case empty ref allele
                yield var_entry.variant_id, i, q
            elif num_pending > 0:
                # allele is not resolved: re-queue
                vqueue.appendleft(var_entry)
                break
            # else: allele does match. discard and continue

    # After last cigar operation, yield ALL resolved variants, pop unresolved variants
    for var_entry in vqueue:
        resolved = list(var_entry.get_resolved())
        num_resolved = len(resolved)
        num_pending = len(var_entry.get_pending())
        if num_resolved >= 1 and num_pending == 0:
            # multiple alleles possible: yield longest
            lengths = [var_entry.alleles[r].length for r in resolved]
            i = resolved[lengths.index(max(lengths))]
            a = var_entry.alleles[i]
            q = a.quality // a.length if a.length > 0 else 30  # Corner case empty ref allele
            yield var_entry.variant_id, i, q

def _detect_alleles_match(variant, entry, bam_read, ref_pos, query_pos, length):
    cdef int query_start = entry.query_start
    cdef int op_start = max(0, entry.query_start - query_pos)
    cdef int ops_consumed
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

def _detect_alleles_insertion(variant, entry, bam_read, ref_pos, query_pos, length):
    cdef int query_start = entry.query_start
    cdef int ops_consumed
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

def _detect_alleles_deletion(variant, entry, bam_read, ref_pos, query_pos, length):
    cdef int ops_consumed = 0
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