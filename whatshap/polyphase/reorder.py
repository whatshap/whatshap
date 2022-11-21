import itertools as it
import logging
from collections import defaultdict
from math import log
from enum import Enum
from copy import deepcopy
from pulp import LpProblem, LpVariable, LpMaximize, LpInteger, COIN_CMD

logger = logging.getLogger(__name__)


class ReorderType(Enum):
    NONE = 0
    MULTI = 1
    COLLAPSED = 2
    COLLAPSED_PRE = 3


class ReorderEvent:
    def __init__(self, position, event_type, best_prob, second_prob, haplotypes):
        self.position = position
        self.event_type = event_type
        ratio = best_prob / (best_prob + second_prob) if best_prob + second_prob > 0 else 0
        self.confidence = ratio
        self.haplotypes = sorted(haplotypes[:])


def run_reordering(
    allele_matrix,
    clustering,
    threads,
    haplotypes,
    block_cut_sensitivity,
    error_rate=0.025,
):
    """
    Main method for the reordering stage of the polyploid phasing algorithm. Input:

    allele_matrix -- The fragment matrix to phase
    clustering -- A list of clusters. Each cluster is a list of read ids, indicating which reads it
                  contains. Every read can only be present in one cluster.
    threads -- The computed cluster paths from the threading stage
    haplotypes -- The computed haplotype sequences (as numeric alleles) from the threading stage
    block_cut_sensitivity -- Policy how conversative the block cuts have to be done. 0 is one
                             phasing block no matter what, 5 is very short blocks

    In this phase, ambiguous sites are reordered by using the most likely arrangement according to
    the reads. Ambiguous sites (resolved in this order) are:

    1. Two (or more) haplotypes being threaded through the same cluster but not having the same
       allele according to precomputed haplotypes. The ordering of these alleles does not matter
       for the threading stage, but multiple of such positions should be ordered to match the input
       reads.
    2. Two (or more) haplotypes being threaded through the same cluster but only a subset of them
       leaves the cluster. Since haplotypes have a history, it matters which haplotype is selected
       to leave, even though the coverage-based threading does not distinguish. Alternatively, two
       (or more) haplotypes on different clusters switch to a new one simultaneously. Here it again
       matters which haplotype goes to which, but all combinations have been equally good in the
       threading stage. These two types of ambiguities are resolved in ascending genome position
       order.
    3. Based on uncertainty from the previous two steps and the selected block cut sensitivity,
       the clustering must eventually be split into blocks.

    """

    # determine snp positions inside clusters
    cwise_snps = find_cluster_snps(threads, haplotypes)

    # phase cluster snps
    phase_cluster_snps(
        threads, haplotypes, cwise_snps, clustering, allele_matrix, error_rate, window_size=32
    )

    # find breakpoints where blocks have to be reordered
    breakpoints = find_breakpoints(threads)

    # compute link log likelihoods between consecutive blocks
    lllh = compute_link_likelihoods(
        threads, haplotypes, breakpoints, clustering, allele_matrix, error_rate
    )
    assert len(lllh) == len(breakpoints)

    # compute optimal permutation of threads for each block
    events = permute_blocks(threads, haplotypes, breakpoints, lllh)

    cut_positions, haploid_cuts = compute_cut_positions(threads, block_cut_sensitivity, events)

    logger.debug(f"Cut positions: {cut_positions}")
    ploidy = len(haplotypes)
    for i in range(ploidy):
        logger.debug(f"Cut positions on phase {i}: {haploid_cuts[i]}")

    return (cut_positions, haploid_cuts, threads, haplotypes)


def find_cluster_snps(threads, haplotypes):
    """
    For each cluster, finds the positions, where it has a multiplicity of at least 2 with at least
    2 different alleles. These positions have to be phased within the clusters.
    """
    cwise_snps = defaultdict(list)
    for pos in range(len(threads)):
        clusters = set()
        alleles = defaultdict(set)
        for i, cid in enumerate(threads[pos]):
            clusters.add(cid)
            alleles[cid].add(haplotypes[i][pos])
        for cid in clusters:
            if len(alleles[cid]) >= 2:
                cwise_snps[cid].append(pos)

    return cwise_snps


def phase_cluster_snps(
    threads, haplotypes, cwise_snps, clustering, allele_matrix, error_rate, window_size
):
    """
    For clusters with multiplicity >= 2 on multiple positions: Bring the alleles of the
    precomputed haplotypes in order by using read information.
    The window size determines how many previous positions are taken into account when resolving
    the next one.
    """
    # sort clusters by position of first haplotype
    collapsed = sorted(
        ((cid, sorted(snps)) for cid, snps in cwise_snps.items()), key=lambda x: x[1][0]
    )

    # iterate cluster-wise
    for cid, snps in collapsed:
        # determine the last window_size many heterozyguous positions on the threads of the first snp
        c_slots = [j for j in range(len(threads[snps[0]])) if threads[snps[0]][j] == cid]
        het_pos = []
        i, w = snps[0] - 1, 0
        while w < window_size and i >= 0:
            if len({haplotypes[s][i] for s in c_slots}) > 1:
                het_pos.append(i)
                w += 1
            i -= 1
        het_pos = het_pos[::-1] + snps
        # w = number of preceding het positions, can be lower than window_size close to chr start

        # use submatrix of relevant reads and positions
        submatrix = allele_matrix.extractSubMatrix(het_pos, clustering[cid], True)

        # reconstruct cluster phasing
        for i, pos in enumerate(snps):
            het_idx = i + w
            c_slots = [j for j in range(len(threads[pos])) if threads[pos][j] == cid]

            configs = []
            # enumerate all permutations of haplotypes at current positions
            for perm in it.permutations(c_slots):
                score = 0
                # for each read: determine number of errors for best fit
                for rid in range(len(submatrix)):
                    if submatrix.getAllele(rid, het_idx) < 0:
                        continue
                    if submatrix.getFirstPos(rid) >= het_idx:
                        continue
                    likelihood = 0.0
                    a_priori = 1 / len(c_slots)
                    for slot in range(len(c_slots)):
                        overlap, errors = 0, 0
                        # skip if read has different allele at current pos than tested config
                        if haplotypes[perm[slot]][pos] != submatrix.getAllele(rid, het_idx):
                            continue
                        # else, compare to window_size many previous positions
                        for (j, a) in submatrix.getRead(rid):
                            if j >= het_idx:
                                break
                            overlap += 1
                            if haplotypes[c_slots[slot]][het_pos[j]] != a:
                                errors += 1
                        likelihood += (
                            a_priori
                            * ((1 - error_rate) ** (overlap - errors))
                            * (error_rate**errors)
                        )
                    score += log(likelihood) if likelihood > 0 else -float("inf")
                configs.append((perm, score))

            configs.sort(key=lambda x: -x[1])
            best_perm = configs[0][0]

            # switch alleles at current position, needed as input for next position(s)
            alleles = [haplotypes[best_perm[j]][pos] for j in range(len(best_perm))]
            for j in range(len(c_slots)):
                haplotypes[c_slots[j]][pos] = alleles[j]
        del submatrix


def find_breakpoints(threads):
    """
    Finds positions p such that between p-1 and p there is an ambiguous switch in the haplotype
    threads. Ambiguous means that two or more threads switch clusters simultaneously or that
    two or more threads were on the same cluster c at position p-1 and one of them is on another
    cluster at position p. Returns a dictionary that for each positions returns the indices of
    affected threads.
    """
    ploidy = len(threads[0])
    breakpoints = dict()

    for i in range(1, len(threads)):
        changed_idx = {j for j in range(ploidy) if threads[i - 1][j] != threads[i][j]}
        affected_clusts = {threads[i - 1][j] for j in changed_idx}
        affected_idx = sorted(j for j in range(ploidy) if threads[i - 1][j] in affected_clusts)

        # ambiguous, if at least two affected
        if len(affected_idx) >= 2:
            breakpoints[i] = affected_idx

    return breakpoints


def compute_link_likelihoods(
    threads, haplotypes, breakpoints, clustering, allele_matrix, error_rate
):
    """
    For each breakpoint and for each pair of threads t1 and t2, computes a log likelihood that the
    left side of t1 is linked to the right side of t2 (left/right = before/after the breakpoint).
    Returns a 2D-list, where the first dimension one entry for each breakpoint and the second
    dimension contains ploidy^2 entries, representing the pairwise linkage likelihoods.
    """
    ploidy = len(threads[0])
    lllh = []
    for pos, affected in breakpoints.items():
        # find the last window_size positions (starting from pos - 1),
        # which are heterozyguous among the haplotype group
        left_pos, right_pos = get_heterozygous_pos_for_haps(haplotypes, affected, pos, 32)
        both_pos = sorted(left_pos + right_pos)

        # use submatrix of relevant reads and positions
        affected_clusts = {threads[pos][h] for h in affected}
        rids = filter(
            lambda r: allele_matrix.getFirstPos(r) < pos <= allele_matrix.getLastPos(r),
            [r for cid in affected_clusts for r in clustering[cid]],
        )
        submatrix = allele_matrix.extractSubMatrix(both_pos, list(rids), True)

        # per read, per haplotype and per side (left/right): compute originating likelihood
        left_llh, right_llh = [], []
        for read in submatrix:
            left_l, right_l = [], []
            for h in range(ploidy):
                l_olp, r_olp, l_err, r_err = 0, 0, 0, 0
                for (j, a) in read:
                    p = both_pos[j]
                    error = 0 if a == haplotypes[h][p] else 1
                    if p < pos:
                        l_olp += 1
                        l_err += error
                    else:
                        r_olp += 1
                        r_err += error
                prob = ((1 - error_rate) ** (l_olp - l_err)) * (error_rate**l_err)
                llh = log(prob) if prob > 0 else -float("inf")
                left_l.append(llh)
                prob = ((1 - error_rate) ** (r_olp - r_err)) * (error_rate**r_err)
                llh = log(prob) if prob > 0 else -float("inf")
                right_l.append(llh)
            left_llh.append(left_l)
            right_llh.append(right_l)

        # per permutation of connections, compute combined likelihood over all reads
        perm_llhs = dict()
        for perm in it.permutations(affected):
            left_h = list(range(ploidy))
            right_h = [perm[affected.index(i)] if i in affected else i for i in range(ploidy)]
            perm_llh = 0
            for i, read in enumerate(submatrix):
                read_llh = -float("inf")
                for left, right in zip(left_h, right_h):
                    read_llh = max(read_llh, left_llh[i][left] + right_llh[i][right])
                perm_llh += read_llh
            perm_llhs[perm] = perm_llh

        # if all permutations infeasible, set a dummy score of 0 to avoid -inf as ILP objective
        if max(perm_llhs) == -float("inf"):
            for perm in perm_llhs:
                perm_llhs[perm] = 0

        lllh.append(perm_llhs)

    return lllh


def get_heterozygous_pos_for_haps(haplotypes, subset, pivot_pos, limit=0):
    """
    For a subset of given haplotypes, returns two lists of positions on which these haplotypes
    contain at least two different alleles. The first list contains positions left to the
    provided pivot position (excluding the pivot itself). The second list contains the positions
    right of the pivot (including the pivot itself if heterozygous). If a limit is provided,
    computes at most so many positions per list.
    """

    left, right = [], []
    j = pivot_pos - 1

    while len(left) < limit and j >= 0:
        if len({haplotypes[h][j] for h in subset}) > 1:
            left.append(j)
        j -= 1
    left = left[::-1]

    j = pivot_pos
    while len(right) < limit and j < len(haplotypes[0]):
        if len({haplotypes[h][j] for h in subset}) > 1:
            right.append(j)
        j += 1
    return left, right


def permute_blocks(threads, haplotypes, breakpoints, lllh):

    ploidy = len(threads[0])
    P = list(range(ploidy))
    B = list(range(len(breakpoints)))
    BE = list(range(len(breakpoints) + 1))
    ext_bp = [0] + list(breakpoints.keys()) + [len(threads)]

    # setup ILP model
    model = LpProblem("PermuteBlocks_p{}_m{}".format(ploidy, len(threads)), LpMaximize)

    # x[b][t][h] = 1 if thread t is placed on haplotype h for block b
    x = [
        [[LpVariable("x_{}_{}_{}".format(b, t, h), 0, 1, LpInteger) for h in P] for t in P]
        for b in BE
    ]

    # for now: fix first permutation
    for t in P:
        model += x[0][t][t] == 1

    # y[b][t1][t2] = 1 if thread t1 is linked to thread t2 over breakpoint b
    y = [
        [[LpVariable("y_{}_{}_{}".format(b, t1, t2), 0, 1, LpInteger) for t2 in P] for t1 in P]
        for b in B
    ]

    # z[b][i] = 1 if i-th combination is used to connect haplotypes over breakpoint b
    z = [
        [LpVariable("z_{}_{}".format(b, i), 0, 1, LpInteger) for i in range(len(lllh[b]))]
        for b in B
    ]

    # for each block there must be one-to-one-relationship between threads and haplotypes
    for i in BE:
        for j in P:
            model += sum([x[i][j][k] for k in P]) == 1
            model += sum([x[i][k][j] for k in P]) == 1

    # iff t1 and t2 both cover haplotype h on blocks b and b+1, set relationship
    for b in B:
        affected = breakpoints[ext_bp[b + 1]]
        for t1 in P:
            for t2 in P:
                if (t1 in affected) != (t2 in affected):
                    # if one of two threads is not affected by the breakpoint, never link these
                    model += y[b][t1][t2] == 0
                elif t1 not in affected:  # and t2 not in affected
                    # unaffected threads only link to themselves
                    if t1 == t2:
                        model += y[b][t1][t2] == 1
                    else:
                        model += y[b][t1][t2] == 0
                for h in P:
                    model += x[b][h][t1] + x[b + 1][h][t2] - 1 <= y[b][t1][t2]
            # every thread must be linked to exactly one other thread over every breakpoint
            model += sum([y[b][t1][t2] for t2 in P]) == 1
            model += sum([y[b][t2][t1] for t2 in P]) == 1

    # ensure correct setting of z-variables, depending on choice of y-variables
    z_weights = dict()
    for b in B:
        left = breakpoints[ext_bp[b + 1]]
        assert left == sorted(left)
        # iterate over all permutations of "left" as given by the keys of lllh for breakpoint b
        for i, right in enumerate(lllh[b].keys()):
            # z[b][i] represents a left-to-right mapping for the haplotypes at b
            z_weights[z[b][i]] = lllh[b][right]
            assert set(left) == set(right)
            # if all y-variables of this permutation are set, force z-variable to 1
            model += z[b][i] >= sum(y[b][l][r] for l, r in zip(left, right)) - len(left) + 1
            # if any of the v-variables is unset, force z-variable to 0
            for l, r in zip(left, right):
                model += z[b][i] <= y[b][l][r]
        # exactly one z-variable must be 1 (should be implied by the above equations anyways)
        model += sum(z[b]) == 1

    # add permutation likelihoods for taken thread links as objective
    model += sum([var * weight for (var, weight) in z_weights.items()])

    # solve model
    solver = COIN_CMD(msg=1)
    model.solve(solver)

    assignments = [[0 for _ in P] for _ in BE]
    for b in BE:
        for t in P:
            for h in P:
                if x[b][t][h].varValue > 0.999:
                    assignments[b][t] = h
                    break
            else:
                assert False

    # shuffle threads and haplotypes according to optimal assignments
    threads_copy = deepcopy(threads)
    haplotypes_copy = deepcopy(haplotypes)
    for i, (s, e) in enumerate(zip(ext_bp[:-1], ext_bp[1:])):
        for p in range(s, e):
            for t in P:
                threads[p][t] = threads_copy[p][assignments[i][t]]
                haplotypes[t][p] = haplotypes_copy[assignments[i][t]][p]

    # TODO: Dummy events
    events = []
    for pos, affected in breakpoints.items():
        events.append(ReorderEvent(pos, ReorderType.COLLAPSED, 0.5, 0.25, affected))

    return events


def compute_cut_positions(path, block_cut_sensitivity, events):
    """
    Takes a threading as input and computes on which positions a cut should be made according the
    cut sensitivity. The levels mean:

    0 -- No cuts at all, even if regions are not connected by any reads
    1 -- Only cut, when regions are not connected by any reads (already done in advance)
    2 -- Only cut, when regions are not connected by a sufficient number of reads (also done)
    3-5 -- Cut on reordering events, depending on confidence values and types.

    The list of cut positions contains the first position of every block. Therefore, position 0 is
    always in the cut list. The second return value is a list of cut positions for every haplotype
    individually.
    """

    cuts = [0]
    hap_cuts = [[0] for _ in range(len(path[0]) if path else 0)]

    if block_cut_sensitivity >= 3:
        if block_cut_sensitivity >= 5:
            pre_thrshld = 0.9
            mul_thrshld = 1.0
            gen_thrshld = 1.0
        elif block_cut_sensitivity == 4:
            pre_thrshld = 0.7
            mul_thrshld = 0.9
            gen_thrshld = 0.9
        else:
            pre_thrshld = 0.5
            mul_thrshld = 0.8
            gen_thrshld = 0.8

        for e in events:
            cut = False
            if e.event_type is ReorderType.COLLAPSED_PRE and e.confidence <= pre_thrshld:
                cut = True
            elif e.event_type is ReorderType.MULTI and e.confidence <= mul_thrshld:
                cut = True
            elif e.event_type is ReorderType.COLLAPSED and e.confidence <= gen_thrshld:
                cut = True
            if cut:
                cuts.append(e.position)
                for h in e.haplotypes:
                    hap_cuts[h].append(e.position)

    return cuts, hap_cuts
