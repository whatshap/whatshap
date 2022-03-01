import itertools as it
import logging
from collections import defaultdict
from .core import HaploThreader

logger = logging.getLogger(__name__)


def run_threading(readset, clustering, ploidy, genotypes, block_cut_sensitivity):
    """
    Main method for the threading stage of the polyploid phasing algorithm. Takes the following input:

    readset -- The fragment matrix to phase
    clustering -- A list of clusters. Each cluster is a list of read ids, indicating which reads it contains. Every read can
                  only be present in one cluster.
    ploidy -- Number of haplotypes to phase
    block_cut_sensitivity -- Policy how conversative the block cuts have to be done. 0 is one phasing block no matter what, 5
                             is very short blocks

    For every variant, the threading algorithm finds a tuple of clusters through which the haplotypes can be threaded with
    minimal cost. Costs arise when the positional coverage of a cluster does not match the number of haplotypes threaded through it
    or when haplotypes switch the cluster on two consecutive positions.
    """

    # compute auxiliary data
    index, rev_index = get_position_map(readset)
    num_vars = len(rev_index)
    positions = get_cluster_start_end_positions(readset, clustering, index)
    coverage = get_coverage(readset, clustering, index)
    cov_map = get_pos_to_clusters_map(coverage, ploidy)
    consensus = get_local_cluster_consensus(readset, clustering, cov_map, positions)

    # compute threading through the clusters
    path = compute_threading_path(
        readset, clustering, num_vars, coverage, cov_map, consensus, ploidy, genotypes
    )

    # we can look at the sequences again to use the most likely continuation, when two or more clusters switch at the same position
    num_clusters = len(clustering)
    c_to_c_global = compute_cluster_to_cluster_similarity(
        readset, clustering, index, consensus, cov_map
    )
    path = improve_path_on_multiswitches(path, num_clusters, c_to_c_global)

    # we can look at the sequences again to use the most likely continuation, when a haplotype leaves a collapsed cluster (currently inactive)
    path = improve_path_on_collapsedswitches(path, num_clusters, c_to_c_global)

    cut_positions, haploid_cuts = compute_cut_positions(path, block_cut_sensitivity, num_clusters)

    logger.debug("Cut positions: {}".format(cut_positions))
    for i in range(ploidy):
        logger.debug("Cut positions on phase {}: {}".format(i, haploid_cuts[i]))

    # compute haplotypes
    haplotypes = []
    for i in range(ploidy):
        hap = ""
        alleles_as_strings = []
        for pos in range(len(path)):
            c_id = path[pos][i]
            allele = consensus[pos][c_id] if c_id in consensus[pos] else -1
            if allele == -1:
                alleles_as_strings.append("n")
            else:
                alleles_as_strings.append(str(allele))
        haplotypes.append(hap.join(alleles_as_strings))

    return (cut_positions, haploid_cuts, path, haplotypes)


def compute_threading_path(
    readset,
    clustering,
    num_vars,
    coverage,
    cov_map,
    consensus,
    ploidy,
    genotypes,
    switch_cost=32.0,
    affine_switch_cost=8.0,
):
    """
    Runs the threading algorithm for the haplotypes using the given costs for switches. The normal switch cost is the
    cost per haplotype, that switches clusters from one position to another (non matching coverage on single position
    is always cost 1.0). The affine switch cost is an additional offset for every position, where a switch occurs.
    These additional costs encourage the threading algorithm to summarize multiple switches on consecutive positions
    into one big switch.
    """

    logger.debug("Computing threading paths ..")

    # arrange data
    compressed_coverage = (
        []
    )  # rearrange the content of "coverage" in a position-wise way, such that the i-th coverage refers to the i-th cluster in cov_map
    compressed_consensus = (
        []
    )  # for every position, give a list of consensi, such that the i-th consensus refers to the i-th cluster in cov_map
    for pos in range(num_vars):
        coverage_list = []
        consensus_list = []
        for i in range(len(cov_map[pos])):
            coverage_list.append(coverage[pos][cov_map[pos][i]])
            consensus_list.append(consensus[pos][cov_map[pos][i]])
        compressed_coverage.append(coverage_list)
        compressed_consensus.append(consensus_list)

    # run threader
    threader = HaploThreader(
        ploidy, switch_cost, affine_switch_cost, True, 16 * 2 ** ploidy if ploidy > 6 else 0
    )
    path = threader.computePathsBlockwise(
        [0], cov_map, compressed_coverage, compressed_consensus, genotypes
    )
    assert len(path) == num_vars

    return path


def compute_cut_positions(path, block_cut_sensitivity, num_clusters):
    """
    Takes a threading as input and computes on which positions a cut should be made according the cut sensitivity. The levels mean:

    0 -- No cuts at all, even if regions are not connected by any reads
    1 -- Only cut, when regions are not connected by any reads (is already done in advance, so nothing to do here)
    2 -- Only cut, when regions are not connected by a sufficient number of reads (also lready done in advance)
    3 -- Cut between two positions, if at least two haplotypes switch their cluster on this transition. In this case it is ambiguous,
         how the haplotype are continued and we might get switch errors if we choose arbitrarily.
    4 -- Additionally to 3, cut every time a haplotype leaves a collapsed region. Collapsed regions are positions, where multiple
         haplotypes go through the same cluster. If one cluster leaves (e.g. due to decline in cluster coverage), we do not know
         which to pick, so we are conservative here. Exception: If a cluster contains a set of multiple haplotypes since the start of
         the current block and hap wants to leave, we do not need a cut, since it does not matter which haps leaves. Default option.
    5 -- Cut every time a haplotype switches clusters. Most conservative, but also very short blocks.

    The list of cut positions contains the first position of every block. Therefore, position 0 is always in the cut list. The second
    return value is a list of cut positions for every haplotype individually.
    """

    cut_positions = [0]
    haploid_cut_positions = []

    if len(path) == 0:
        return cut_positions

    ploidy = len(path[0])
    haploid_cut_positions = [[0] for _ in range(ploidy)]

    dissim_threshold = 1
    rise_fall_dissim = 0
    if block_cut_sensitivity >= 3:
        if block_cut_sensitivity >= 5:
            # cut every time a haplotype jumps
            dissim_threshold = 1
            rise_fall_dissim = ploidy + 1
        elif block_cut_sensitivity == 4:
            # cut for every multi-switch and for every rise-fall-ploidy change
            dissim_threshold = 2
            rise_fall_dissim = ploidy + 1
        else:
            # cut for every multi-jump
            dissim_threshold = 2
            rise_fall_dissim = 0

    if block_cut_sensitivity >= 3:
        copynrs = []
        for i in range(0, len(path)):
            copynr = defaultdict(int)
            for j in range(0, ploidy):
                if path[i][j] not in copynr:
                    copynr[path[i][j]] = 0
                copynr[path[i][j]] += 1
            copynrs.append(copynr)

        cpn_rising = [False for c_id in range(num_clusters)]

        for i in range(1, len(path)):
            dissim = 0
            clusters_cut = set()
            for j in range(0, ploidy):
                old_c = path[i - 1][j]
                new_c = path[i][j]
                if old_c != new_c:
                    clusters_cut.add(old_c)
                    rise_fall = False
                    # check if previous cluster went down from copy number >= 2 to a smaller one >= 1
                    if copynrs[i - 1][old_c] > copynrs[i][old_c] >= 1:
                        if cpn_rising[old_c]:
                            rise_fall = True
                    # check if new cluster went up from copy number >= 1 to a greater one >= 2
                    if copynrs[i][new_c] > copynrs[i - 1][new_c] >= 1:
                        cpn_rising[new_c] = True
                    # check if one cluster has been rising and then falling in the current block
                    if rise_fall:
                        dissim += rise_fall_dissim

                    # count general switches
                    dissim += 1

            # TODO: Avoid deep indentation by using functions
            if dissim >= dissim_threshold:
                cpn_rising = [False] * num_clusters
                cut_positions.append(i)

                # get all cut threads
                threads_cut = [j for j in range(ploidy) if path[i - 1][j] in clusters_cut]
                for thread in threads_cut:
                    haploid_cut_positions[thread].append(i)

    return cut_positions, haploid_cut_positions


def compute_cluster_to_cluster_similarity(readset, clustering, index, consensus, cov_map):
    """
    For every position p, compute the similarity between present clusters at position p-1 and
    clusters at position p. Format is: c_to_c_sim[position][(cluster_p-1, cluster_p)]
    """
    num_vars = len(consensus)
    num_clusters = len(clustering)
    coverage_abs = get_coverage_absolute(readset, clustering, index)
    c_to_c_sim = [defaultdict(float) for _ in range(num_vars)]

    cluster_zeroes = [dict() for c_id in range(num_clusters)]
    cluster_ones = [dict() for c_id in range(num_clusters)]
    for pos in range(num_vars):
        for c_id in consensus[pos]:
            cluster_zeroes[c_id][pos] = coverage_abs[pos][c_id] * (1 - consensus[pos][c_id])
            cluster_ones[c_id][pos] = coverage_abs[pos][c_id] * consensus[pos][c_id]

    for var in range(1, num_vars):
        for i, c1 in enumerate(cov_map[var - 1]):
            for j, c2 in enumerate(cov_map[var]):
                same = 0
                diff = 0
                # Use a sliding window of positions as basis for consensus similarity
                for pos in range(max(0, var - 10), min(num_vars - 1, var + 9)):
                    if pos in cluster_zeroes[c1] and pos in cluster_zeroes[c2]:
                        same += (
                            cluster_zeroes[c1][pos] * cluster_zeroes[c2][pos]
                            + cluster_ones[c1][pos] * cluster_ones[c2][pos]
                        )
                        diff += (
                            cluster_zeroes[c1][pos] * cluster_ones[c2][pos]
                            + cluster_ones[c1][pos] * cluster_zeroes[c2][pos]
                        )
                c_to_c_sim[var][(c1, c2)] = same / (same + diff) if same > 0 else 0

    return c_to_c_sim


def improve_path_on_multiswitches(path, num_clusters, cluster_sim):
    """
    Post processing step after the threading. If two or more haplotypes switch clusters on the same position, we could use
    the similarity scores between the clusters to find the most likely continuation of the switching haplotypes. See the
    description of the compute_cut_positions method for more details about block cuts.
    """

    if len(path) == 0:
        return []

    corrected_path = []
    corrected_path.append(path[0])
    ploidy = len(path[0])
    current_perm = tuple(range(ploidy))
    invers_perm = [i for i in range(ploidy)]

    for i in range(1, len(path)):
        changed = []  # set of haplotypes, that changed cluster at current position
        for j in range(0, ploidy):
            old_c = path[i - 1][j]
            new_c = path[i][j]
            if old_c != new_c:
                # count general switches
                changed.append(j)

        if len(changed) >= 2:
            # if at least two threads changed cluster: find optimal permutation of changed clusters
            left_c = [path[i - 1][j] for j in changed]
            right_c = [path[i][j] for j in changed]
            actual_score = sum(
                [cluster_sim[i][(left_c[j], right_c[j])] for j in range(len(changed))]
            )
            best_score = actual_score
            best_perm = tuple(range(len(changed)))
            for perm in it.permutations(range(len(changed))):
                score = 0
                for j, left in enumerate(left_c):
                    score += cluster_sim[i][(left, right_c[perm[j]])]
                if score > best_score:
                    best_score = score
                    best_perm = perm

            # apply local best permutation to current global permutation
            current_perm_copy = list(current_perm)
            for j in range(len(changed)):
                current_perm_copy[changed[j]] = current_perm[changed[best_perm[j]]]
            current_perm = tuple(current_perm_copy)
            for j in range(ploidy):
                invers_perm[current_perm[j]] = j

        # apply current optimal permutation to local cluster config and add to corrected path
        corrected_path.append([path[i][j] for j in invers_perm])

    return corrected_path


def improve_path_on_collapsedswitches(path, num_clusters, cluster_sim):
    """
    Post processing step after the threading. If a haplotype leaves a cluster on a collapsed region, we could use the
    similarity scores between the clusters to find the most likely continuation of the switching haplotypes. See the
    description of the compute_cut_positions method for more details about block cuts.
    """
    if len(path) == 0:
        return []

    corrected_path = []
    corrected_path.append(path[0])
    ploidy = len(path[0])
    current_perm = tuple(range(ploidy))
    invers_perm = [i for i in range(ploidy)]

    copynrs = []
    for i in range(0, len(path)):
        copynr = defaultdict(int)
        for j in range(0, ploidy):
            if path[i][j] not in copynr:
                copynr[path[i][j]] = 0
            copynr[path[i][j]] += 1
        copynrs.append(copynr)

    for i in range(1, len(path)):
        changed = []
        # iterate over present cluster ids
        for c_id in copynrs[i]:
            if copynrs[i - 1][c_id] >= 2:
                # for all collapsed clusters: find haplotypes, which go through and check whether one of them exits
                outgoing_c = False
                affected = []
                for j in range(ploidy):
                    if path[i - 1][j] == c_id:
                        affected.append(j)
                        if path[i][j] != c_id:
                            outgoing_c = True
                # if haplotypes leaves collapsed cluster, all other might be equally suited, so add them
                if outgoing_c:
                    changed.append(affected)

        for h_group in changed:
            # for every group of haplotypes coming from a collapsed cluster:
            collapsed_cid = path[i - 1][h_group[0]]
            left_c = []

            # find last cluster before collapsed one for every haplotype (or use collapsed one if this does not exist)
            for j in h_group:
                pos = i - 1
                while pos >= 0:
                    if path[pos][j] != collapsed_cid:
                        left_c.append(path[pos][j])
                        break
                    else:
                        pos -= 1
                if pos == -1:
                    left_c.append(collapsed_cid)
            right_c = [path[i][j] for j in h_group]

            # we need to catch the case, where we compare a cluster with itself
            ident_sim = 0
            for c1 in left_c:
                for c2 in right_c:
                    if c1 != c2:
                        ident_sim = max(ident_sim, cluster_sim[i][(c1, c2)])
            ident_sim = ident_sim * 2 + 1

            for j in range(len(h_group)):
                actual_score = sum(
                    [
                        cluster_sim[i][(left_c[j], right_c[j])]
                        if left_c[j] != right_c[j]
                        else ident_sim
                        for j in range(len(h_group))
                    ]
                )
            best_score = actual_score
            best_perm = tuple(range(len(h_group)))
            for perm in it.permutations(range(len(h_group))):
                score = 0
                for j, left in enumerate(left_c):
                    score += (
                        cluster_sim[i][(left, right_c[perm[j]])]
                        if left != right_c[perm[j]]
                        else ident_sim
                    )
                if score > best_score:
                    best_score = score
                    best_perm = perm

            # apply local best permutation to current global permutation
            current_perm_copy = list(current_perm)
            for j in range(len(h_group)):
                current_perm_copy[h_group[j]] = current_perm[h_group[best_perm[j]]]
            current_perm = tuple(current_perm_copy)
            for j in range(ploidy):
                invers_perm[current_perm[j]] = j

        # apply current optimal permutation to local cluster config and add to corrected path
        corrected_path.append([path[i][j] for j in invers_perm])

    return corrected_path


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


def get_coverage(readset, clustering, pos_index):
    """
    Returns a list, which for every position contains a dictionary, mapping a cluster id to
    a relative coverage on this position.
    """
    num_vars = len(pos_index)
    num_clusters = len(clustering)
    coverage = [dict() for pos in range(num_vars)]
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


def get_coverage_absolute(readset, clustering, pos_index):
    """
    Returns a list, which for every position contains a dictionary, mapping a cluster id to
    an absolute coverage on this position.
    """
    num_vars = len(pos_index)
    num_clusters = len(clustering)
    coverage = [dict() for pos in range(num_vars)]
    for c_id in range(num_clusters):
        for read in clustering[c_id]:
            for pos in [pos_index[var.position] for var in readset[read]]:
                if c_id not in coverage[pos]:
                    coverage[pos][c_id] = 0
                coverage[pos][c_id] += 1

    return coverage


def get_cluster_start_end_positions(readset, clustering, pos_index):
    num_clusters = len(clustering)
    positions = {}
    for c_id in range(num_clusters):
        read = clustering[c_id][0]
        start = pos_index[readset[read][0].position]
        end = pos_index[readset[read][-1].position]
        for read in clustering[c_id]:
            readstart = pos_index[readset[read][0].position]
            readend = pos_index[readset[read][-1].position]
            if readstart < start:
                start = readstart
            if readend > end:
                end = readend
        positions[c_id] = (start, end)
    assert len(positions) == num_clusters
    return positions


def get_local_cluster_consensus(readset, clustering, cov_map, positions):
    """
    Returns a list, which for every position contains a dictionary, mapping a cluster id to
    its consensus on this position.
    """
    return [
        {c_id: pos_cons[c_id][0] for c_id in pos_cons}
        for pos_cons in get_local_cluster_consensus_withfrac(
            readset, clustering, cov_map, positions
        )
    ]


def get_local_cluster_consensus_withfrac(readset, clustering, cov_map, positions):
    # Map genome positions to [0,l)
    index = {}
    rev_index = []
    num_vars = 0
    for position in readset.get_positions():
        index[position] = num_vars
        rev_index.append(position)
        num_vars += 1

    relevant_pos = [[] for i in range(len(clustering))]
    for pos in range(num_vars):
        for c in cov_map[pos]:
            relevant_pos[c].append(pos)

    clusterwise_consensus = [
        get_single_cluster_consensus_frac(readset, clustering[i], index, relevant_pos[i])
        for i in range(len(clustering))
    ]
    whole_consensus = []
    for pos in range(num_vars):
        newdict = defaultdict()
        for c in cov_map[pos]:
            newdict[c] = clusterwise_consensus[c][pos]
        whole_consensus.append(newdict)
    return whole_consensus


def get_single_cluster_consensus_frac(readset, cluster, index, relevant_pos):
    # Count zeroes and one for every position
    poswise_allelecount = dict()
    for read in cluster:
        for var in readset[read]:
            pos = index[var.position]
            if pos not in poswise_allelecount:
                poswise_allelecount[pos] = dict()
            if var.allele not in poswise_allelecount[pos]:
                poswise_allelecount[pos][var.allele] = 0
            poswise_allelecount[pos][var.allele] += 1

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
