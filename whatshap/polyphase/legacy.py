import logging
from collections import defaultdict

from whatshap.polyphase.solver import AlleleMatrix

logger = logging.getLogger(__name__)

# Scipy renamed the binom_test method for version >=1.12. We avoid requirement bumping for this
# legacy module
try:
    from scipy.stats import binom_test
except ImportError:
    from scipy.stats import binomtest


def find_inconsistencies(am: AlleleMatrix, clustering, ploidy):
    # Returns the number of cluster positions with inconsistencies
    # (counts position multiple times, if multiple clusters are inconsistent there)
    # Also returns a list of read pairs, which need to be seperated
    num_inconsistent_positions = 0
    separated_pairs = []
    exp_error = 0.05
    p_val_threshold = 0.02

    # Compute consensus and coverage
    index, rev_index = get_position_map(am)
    num_vars = am.getNumPositions()

    coverage = get_coverage(am, clustering)
    cov_map = get_pos_to_clusters_map(coverage, ploidy)
    positions = get_cluster_start_end_positions(am, clustering)
    abs_coverage = get_coverage_absolute(am, clustering)
    consensus = get_local_cluster_consensus_withfrac(am, clustering, cov_map, positions)

    # Search for positions in clusters with ambivalent consensus
    for pos in range(num_vars):
        for c_id in coverage[pos]:
            if c_id not in consensus[pos]:
                continue
            # do binomial hypothesis test, whether the deviations from majority allele is significant enough for splitting
            abs_count = abs_coverage[pos][c_id]
            abs_dev = int(abs_count * (1 - consensus[pos][c_id][1]))
            # ugly, but found no more elegant way yet to detect which import was successful
            try:
                p_val = binom_test(abs_dev, abs_count, exp_error, alternative="greater")
            except NameError:
                p_val = binomtest(abs_dev, abs_count, exp_error, alternative="greater").pvalue
            if p_val < p_val_threshold:
                num_inconsistent_positions += 1
                zero_reads = []
                one_reads = []
                for read in clustering[c_id]:
                    for rpos, allele in am.getRead(read):
                        if rpos == pos:
                            if allele == 0:
                                zero_reads.append(read)
                            else:
                                one_reads.append(read)
                for r0 in zero_reads:
                    for r1 in one_reads:
                        separated_pairs.append((r0, r1))

    return num_inconsistent_positions, separated_pairs


def get_position_map(am):
    """
    Returns a mapping of genome (bp) positions to virtual positions (from 0 to l).
    """
    # Map genome positions to [0,l)
    index = {}
    rev_index = []
    num_vars = 0

    for position in am.getPositions():
        index[position] = num_vars
        rev_index.append(position)
        num_vars += 1

    return index, rev_index


def get_pos_to_clusters_map(coverage, ploidy):
    """
    For every position, computes a list of relevant clusters for the threading
    algorithm. Relevant means, that the relative coverage is at least 1/8 of
    what a single haplotype is expected to have for the given ploidy. Apart
    from that, at least <ploidy> and at most <2*ploidy> many clusters are
    selected to avoid exponential blow-up.
    """
    cov_map = [[] for _ in range(len(coverage))]
    for pos in range(len(coverage)):
        sorted_cids = sorted(
            [cid for cid in coverage[pos]], key=lambda x: coverage[pos][x], reverse=True
        )
        cut_off = min(len(sorted_cids), 2 * ploidy)
        for i in range(ploidy, min(len(sorted_cids), 2 * ploidy)):
            if coverage[pos][sorted_cids[i]] < (1.0 / (8.0 * ploidy)):
                cut_off = i
                break
        cov_map[pos] = sorted_cids[:cut_off]
    return cov_map


def get_coverage(am, clustering):
    """
    Returns a list, which for every position contains a dictionary, mapping a cluster id to
    a relative coverage on this position.
    """
    num_vars = am.getNumPositions()
    num_clusters = len(clustering)
    coverage = [dict() for pos in range(num_vars)]
    coverage_sum = [0 for pos in range(num_vars)]
    for c_id in range(num_clusters):
        for read in clustering[c_id]:
            for pos, allele in am.getRead(read):
                if c_id not in coverage[pos]:
                    coverage[pos][c_id] = 0
                coverage[pos][c_id] += 1
                coverage_sum[pos] += 1

    for pos in range(num_vars):
        for c_id in coverage[pos]:
            coverage[pos][c_id] = coverage[pos][c_id] / coverage_sum[pos]

    return coverage


def get_coverage_absolute(am, clustering):
    """
    Returns a list, which for every position contains a dictionary, mapping a cluster id to
    an absolute coverage on this position.
    """
    num_vars = am.getNumPositions()
    num_clusters = len(clustering)
    coverage = [dict() for pos in range(num_vars)]
    for c_id in range(num_clusters):
        for read in clustering[c_id]:
            for pos, allele in am.getRead(read):
                if c_id not in coverage[pos]:
                    coverage[pos][c_id] = 0
                coverage[pos][c_id] += 1

    return coverage


def get_cluster_start_end_positions(am, clustering):
    num_clusters = len(clustering)
    positions = {}
    for c_id in range(num_clusters):
        read = clustering[c_id][0]
        start = am.getFirstPos(read)
        end = am.getLastPos(read)
        for read in clustering[c_id]:
            start = min(start, am.getFirstPos(read))
            end = max(end, am.getLastPos(read))
        positions[c_id] = (start, end)
    assert len(positions) == num_clusters
    return positions


def get_local_cluster_consensus_withfrac(am, clustering, cov_map, positions):
    # Map genome positions to [0,l)
    num_vars = am.getNumPositions()
    relevant_pos = [[] for i in range(len(clustering))]
    for pos in range(num_vars):
        for c in cov_map[pos]:
            relevant_pos[c].append(pos)

    clusterwise_consensus = [
        get_single_cluster_consensus_frac(am, clustering[i], relevant_pos[i])
        for i in range(len(clustering))
    ]
    whole_consensus = []
    for pos in range(num_vars):
        newdict = defaultdict()
        for c in cov_map[pos]:
            newdict[c] = clusterwise_consensus[c][pos]
        whole_consensus.append(newdict)
    return whole_consensus


def get_single_cluster_consensus_frac(am, cluster, relevant_pos):
    # Count zeroes and one for every position
    poswise_allelecount = dict()
    for read in cluster:
        for pos, allele in am.getRead(read):
            if pos not in poswise_allelecount:
                poswise_allelecount[pos] = dict()
            if allele not in poswise_allelecount[pos]:
                poswise_allelecount[pos][allele] = 0
            poswise_allelecount[pos][allele] += 1

    # Determine majority allele
    cluster_consensus = {}
    for pos in relevant_pos:
        if pos in poswise_allelecount:
            max_allele = 0
            max_count = 0
            sum_count = 0
            for allele in sorted(poswise_allelecount[pos]):
                cur_count = poswise_allelecount[pos][allele]
                sum_count += cur_count
                if cur_count > max_count:
                    max_allele = allele
                    max_count = cur_count
            cluster_consensus[pos] = (max_allele, max_count / sum_count)
        else:
            cluster_consensus[pos] = (0, 1.0)

    return cluster_consensus
