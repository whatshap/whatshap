import logging
from collections import defaultdict
from dataclasses import dataclass
from queue import Queue
from typing import List

from pulp import listSolvers, getSolver
from whatshap.core import ReadSet
from whatshap.polyphase.solver import AlleleMatrix

logger = logging.getLogger(__name__)


class SolverError(Exception):
    pass


@dataclass
class PolyphaseParameter:
    ploidy: int
    ce_bundle_edges: bool
    distrust_genotypes: bool
    min_overlap: int
    block_cut_sensitivity: int
    plot_clusters: bool
    plot_threading: bool
    threads: int
    use_prephasing: bool


class PhaseBreakpoint:
    def __init__(self, position, haplotypes, confidence):
        self.position = position
        self.haplotypes = sorted(haplotypes[:])
        self.confidence = confidence


@dataclass
class PolyphaseBlockResult:
    block_id: int
    clustering: List[List[int]]
    threads: List[List[int]]
    haplotypes: List[int]
    breakpoints: List[PhaseBreakpoint]


@dataclass
class PolyphaseResult:
    clustering: List[List[int]]
    threads: List[List[int]]
    haplotypes: List[int]
    breakpoints: List[PhaseBreakpoint]


def get_coverage(allele_matrix, clustering):
    """
    Returns a list, which for every position contains a dictionary, mapping a cluster id to
    a relative coverage on this position.
    """
    num_vars = allele_matrix.getNumPositions()
    num_clusters = len(clustering)
    coverage = [defaultdict(int) for pos in range(num_vars)]
    coverage_sum = [0 for pos in range(num_vars)]
    for c_id in range(num_clusters):
        for read in clustering[c_id]:
            for pos, allele in allele_matrix.getRead(read):
                if c_id not in coverage[pos]:
                    coverage[pos][c_id] = 0
                coverage[pos][c_id] += 1
                coverage_sum[pos] += 1

    for pos in range(num_vars):
        for c_id in coverage[pos]:
            coverage[pos][c_id] = coverage[pos][c_id] / coverage_sum[pos]

    return coverage


def compute_block_starts(am, ploidy, single_linkage=False):
    """
    Based on the connectivity of the reads, we want to divide the phasing input, as non- or poorly
    connected regions can be phased independently. This is done based on how pairs of variants are
    connected. There are two modes how to decide whether two variants are connected:

    single_linkage=True -- If there exists a read in the allele matrix which covers both variants,
                           they count as connected
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

    num_vars = am.getNumPositions()

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
    logger.debug(f"Cut position threshold: coverage >= {cut_threshold}")

    # start by looking at neighbouring
    link_to_next = [0 for i in range(num_vars)]
    for read in am:
        pos_list = [pos for (pos, allele) in read]
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
    link_coverage = [defaultdict(int) for i in range(num_clust)]
    for i, read in enumerate(am):
        covered_pos_clusts = {pos_clust[pos] for (pos, allele) in read}
        for p1 in covered_pos_clusts:
            for p2 in covered_pos_clusts:
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


def create_genotype_list(variant_table, sample):
    """
    Creates a list, which stores a dictionary for every position. The dictionary maps every allele
    to its frequency in the genotype of the respective position.
    """
    all_genotypes = variant_table.genotypes_of(sample)
    genotype_list = []
    for pos in range(len(all_genotypes)):
        allele_count = dict()
        for allele in all_genotypes[pos].as_vector():
            if allele not in allele_count:
                allele_count[allele] = 0
            allele_count[allele] += 1
        genotype_list.append(allele_count)
    return genotype_list


def extract_partial_phasing(variant_table, sample, ploidy):
    readset = ReadSet()
    vars = variant_table.variants
    for read in variant_table.phased_blocks_as_reads(sample, vars, 0, 0, target_ploidy=ploidy):
        readset.add(read)
    if len(readset) > 0:
        am = AlleleMatrix(readset)
        assert len(am) % ploidy == 0
        for i in range(0, len(am), ploidy):
            assert all([am.getFirstPos(i) == am.getFirstPos(i + j) for j in range(1, ploidy)])
            assert all([am.getLastPos(i) == am.getLastPos(i + j) for j in range(1, ploidy)])
        return am
    else:
        return None


def get_ilp_solver():
    """
    Sets up a solver object with suppressed cmd output from available solvers.
    """
    solvers = listSolvers(onlyAvailable=True)
    prio = ["GUROBI_CMD", "GUROBI", "COIN_CMD", "PULP_CBC_CMD"]
    for name in prio:
        if name in solvers:
            return getSolver(name, msg=0)
    if len(solvers) > 0:
        return getSolver(solvers[0], msg=0)
    else:
        raise SolverError("No ILP solver is available for PuLP.")
