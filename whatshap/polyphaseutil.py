import logging

from collections import defaultdict
from .core import Read, ReadSet
from queue import Queue

logger = logging.getLogger(__name__)


def get_position_map(readset):
    """
    Returns a mapping of genome (bp) positions to virtual positions (from 0 to l).
    """
    # Map genome positions to [0,l)
    index = {}
    rev_index = []
    num_vars = 0

    for position in readset.get_positions():
        index[position] = num_vars
        rev_index.append(position)
        num_vars += 1

    return index, rev_index


def get_coverage(readset, clustering):
    """
    Returns a list, which for every position contains a dictionary, mapping a cluster id to
    a relative coverage on this position.
    """
    pos_index, rev_index = get_position_map(readset)
    num_vars = len(pos_index)
    num_clusters = len(clustering)
    coverage = [defaultdict(int) for pos in range(num_vars)]
    coverage_sum = [0 for pos in range(num_vars)]
    for c_id in range(num_clusters):
        for read in clustering[c_id]:
            for pos in [pos_index[var.position] for var in readset[read]]:
                if c_id not in coverage[pos]:
                    coverage[pos][c_id] = 0
                coverage[pos][c_id] += 1
                coverage_sum[pos] += 1

    for pos in range(num_vars):
        for c_id in coverage[pos]:
            coverage[pos][c_id] = coverage[pos][c_id] / coverage_sum[pos]

    return coverage


def haplotypes_to_strings(haplotypes):
    haplostrings = []
    for i, h in enumerate(haplotypes):
        haplostrings.append("".join(map(lambda x: str(x) if x >= 0 else "n", h)))
    return haplostrings


def split_readset(readset, ext_block_starts):
    """
    Creates one sub-readset for every block. Reads which cross block borders are also split, so
    parts of them appear in multiple blocks. Reads inside a sub-readset are trimmed, such that they
    do not contain variants outside of their associated blocks.
    """

    index, rev_index = get_position_map(readset)
    var_to_block = [0 for _ in range(ext_block_starts[-1])]
    for i in range(len(ext_block_starts) - 1):
        for var in range(ext_block_starts[i], ext_block_starts[i + 1]):
            var_to_block[var] = i

    block_readsets = [ReadSet() for i in range(len(ext_block_starts) - 1)]
    for i, read in enumerate(readset):
        if not read.is_sorted():
            read.sort()
        start = var_to_block[index[read[0].position]]
        end = var_to_block[index[read[-1].position]]
        if start == end:
            # if read lies entirely in one block, copy it into according readset
            block_readsets[start].add(read)
        else:
            # split read by creating one new read for each covered block
            current_block = start
            read_slice = Read(
                name=read.name,
                source_id=read.source_id,
                sample_id=read.sample_id,
                reference_start=read.sample_id,
                BX_tag=read.BX_tag,
            )
            for variant in read:
                if var_to_block[index[variant.position]] != current_block:
                    block_readsets[current_block].add(read_slice)
                    current_block = var_to_block[index[variant.position]]
                    read_slice = Read(
                        name=str(current_block) + "_" + read.name,
                        source_id=read.source_id,
                        sample_id=read.sample_id,
                        reference_start=read.sample_id,
                        BX_tag=read.BX_tag,
                    )
                read_slice.add_variant(variant.position, variant.allele, variant.quality)
            block_readsets[current_block].add(read_slice)
    return block_readsets


def compute_block_starts(readset, ploidy, single_linkage=False):
    """
    Based on the connectivity of the reads, we want to divide the phasing input, as non- or poorly
    connected regions can be phased independently. This is done based on how pairs of variants are
    connected. There are two modes how to decide whether two variants are connected:

    single_linkage=True -- If there exists a read in the readset, which covers both variants, they
                           are connected
    single_linkage=False -- In order to connect two variants, we need at least reads from ploidy-1
                            different haplotypes. Two variants count as connected, if there
                            sufficiently many reads covering both variants, with "sufficient"
                            meaning, that the connecting reads have a chance of at least 98% that
                            they cover at least ploidy-1 haplotypes.

    First, only consecutive pairs are inspected. Then, this connectivity is made transitive, i.e.
    if the pair (A,C) is connected, as well as the pair (B,C), then (A,B) is also connected. If the
    special case occurs, that for three variants (in this order) A and C are connected, but neither
    is connected to variant B in between them, then the variants are still divided as A|B|C, even
    though A and C are actually connected. This is because the following steps require the splits
    to be intervals with no "holes" inside them.
    """

    pos_index, rev_index = get_position_map(readset)
    num_vars = len(pos_index)

    # special case
    if num_vars == 0:
        return []

    # cut threshold
    if ploidy == 2 or single_linkage:
        cut_threshold = 1
    else:
        cut_threshold = ploidy * ploidy
        for i in range(ploidy - 1, ploidy * ploidy):
            # chance to cover at most ploidy-2 haplotypes with i reads
            cut_threshold = i
            if ploidy * pow((ploidy - 2) / ploidy, i) < 0.02:
                cut_threshold = i
                break
    logger.debug("Cut position threshold: coverage >= {}".format(cut_threshold))

    # start by looking at neighbouring
    link_to_next = [0 for i in range(num_vars)]
    starts = []
    ends = []
    for read in readset:
        pos_list = [pos_index[var.position] for var in read]
        starts.append(pos_list[0])
        ends.append(pos_list[-1])
        for i in range(len(pos_list) - 1):
            if pos_list[i] + 1 == pos_list[i + 1]:
                link_to_next[pos_list[i]] += 1

    pos_clust = [0 for i in range(num_vars)]
    for i in range(1, num_vars):
        if link_to_next[i - 1] >= cut_threshold:
            pos_clust[i] = pos_clust[i - 1]
        else:
            pos_clust[i] = pos_clust[i - 1] + 1
    num_clust = pos_clust[-1] + 1

    # find linkage between clusters
    link_coverage = [dict() for i in range(num_clust)]
    for i, (start, end) in enumerate(zip(starts, ends)):
        covered_pos_clusts = set([pos_clust[pos_index[var.position]] for var in readset[i]])
        for p1 in covered_pos_clusts:
            for p2 in covered_pos_clusts:
                if p2 not in link_coverage[p1]:
                    link_coverage[p1][p2] = 0
                link_coverage[p1][p2] += 1

    # merge clusters
    merged_clust = [-1 for i in range(num_clust)]
    new_num_clust = 0
    for i in range(num_clust):
        if merged_clust[i] >= 0:
            continue
        q = Queue()
        q.put(i)
        merged_clust[i] = new_num_clust
        while not q.empty():
            cur = q.get()
            for linked in link_coverage[cur]:
                if merged_clust[linked] < 0 and link_coverage[cur][linked] >= cut_threshold:
                    q.put(linked)
                    merged_clust[linked] = new_num_clust
        new_num_clust += 1

    # determine cut positions
    cuts = [0]
    for i in range(1, num_vars):
        if merged_clust[pos_clust[i]] != merged_clust[pos_clust[i - 1]]:
            cuts.append(i)

    return cuts


def create_genotype_list(phasable_variant_table, sample):
    """
    Creates a list, which stores a dictionary for every position. The dictionary maps every allele
    to its frequency in the genotype of the respective position.
    """
    all_genotypes = phasable_variant_table.genotypes_of(sample)
    genotype_list = []
    for pos in range(len(all_genotypes)):
        allele_count = dict()
        for allele in all_genotypes[pos].as_vector():
            if allele not in allele_count:
                allele_count[allele] = 0
            allele_count[allele] += 1
        genotype_list.append(allele_count)
    return genotype_list
