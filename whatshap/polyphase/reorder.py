import itertools as it
import logging
from collections import defaultdict
from bisect import bisect_left
from math import log, exp
from enum import Enum
from copy import deepcopy
from pulp import LpProblem, LpVariable, LpMaximize, LpInteger, PULP_CBC_CMD

logger = logging.getLogger(__name__)


class ReorderType(Enum):
    NONE = 0
    MULTI = 1
    COLLAPSED = 2
    COLLAPSED_PRE = 3


class PhaseBreakpoint:
    def __init__(self, position, haplotypes, confidence):
        self.position = position
        self.haplotypes = sorted(haplotypes[:])
        self.confidence = confidence


def run_reordering(
    allele_matrix,
    clustering,
    threads,
    haplotypes,
    prephasing,
    block_cut_sensitivity,
    error_rate=0.05,
):
    """
    Main method for the reordering stage of the polyploid phasing algorithm. Input:

    allele_matrix -- The fragment matrix to phase
    clustering -- A list of clusters. Each cluster is a list of read ids, indicating which reads it
                  contains. Every read can only be present in one cluster.
    threads -- The computed cluster paths from the threading stage
    haplotypes -- The computed haplotype sequences (as numeric alleles) from the threading stage
    prephasing -- Positions for which a phasing already exists, encoded as reads
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
       threading stage. These two types of ambiguities are resolved via an ILP model based on read
       supporting for every possible configuration.
    3. Based on uncertainty from the previous two steps and the selected block cut sensitivity,
       the clustering must eventually be split into blocks.

    """
    ploidy = len(haplotypes)
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

    if prephasing:
        aff = compute_phase_affiliation(
            allele_matrix, haplotypes, breakpoints, prephasing, error_rate
        )
    else:
        aff = None

    # compute optimal permutation of threads for each block
    perms = get_optimal_permutations(breakpoints, lllh, ploidy, aff)
    permute_blocks(threads, haplotypes, breakpoints, lllh, perms)

    return threads, haplotypes, breakpoints


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
                    a_priori = 1 / len(c_slots)
                    # count overlaps and errors first
                    hits, errors = [], []
                    for slot in range(len(c_slots)):
                        overlap, error = 0, 0
                        # skip if read has different allele at current pos than tested config
                        if haplotypes[perm[slot]][pos] != submatrix.getAllele(rid, het_idx):
                            continue
                        # else, compare to window_size many previous positions
                        for (j, a) in submatrix.getRead(rid):
                            if j >= het_idx:
                                break
                            overlap += 1
                            if haplotypes[c_slots[slot]][het_pos[j]] != a:
                                error += 1
                        hits.append(overlap - error)
                        errors.append(error)
                    # use log-space as early as possible to avoid float underflow
                    err_pivot = min(errors)
                    hit_pivot = hits[errors.index(err_pivot)]
                    likelihood = 0.0
                    for hit, error in zip(hits, errors):
                        likelihood += (
                            a_priori
                            * ((1 - error_rate) ** (hit - hit_pivot))
                            * (error_rate ** (error - err_pivot))
                        )
                    """
                    log(w_1 * x^(y_1) + w_2 * x^(y_2)) =
                    log(w_1 * x^(y_1 - y) + w_2 * x^(y_2 - y)) + y * log(x)
                    for w_1 + w_2 = 1
                    """
                    score += log(likelihood) if likelihood > 0 else -float("inf")
                    score += hit_pivot * log(1 - error_rate) + err_pivot * log(error_rate)
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
    breakpoints = []

    for i in range(1, len(threads)):
        changed_idx = {j for j in range(ploidy) if threads[i - 1][j] != threads[i][j]}
        affected_clusts = {threads[i - 1][j] for j in changed_idx}
        affected_haps = sorted(j for j in range(ploidy) if threads[i - 1][j] in affected_clusts)

        # ambiguous, if at least two affected
        if len(affected_haps) >= 2:
            breakpoints.append(PhaseBreakpoint(i, affected_haps, 0.0))

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
    for pos, affected in [(b.position, b.haplotypes) for b in breakpoints]:  # breakpoints.items():
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
                llh = log(1 - error_rate) * (l_olp - l_err) + log(error_rate) * l_err
                left_l.append(llh)
                llh = log(1 - error_rate) * (r_olp - r_err) + log(error_rate) * r_err
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

        assert max(perm_llhs.values()) > -float("inf")
        lllh.append(perm_llhs)

    assert len(lllh) == len(breakpoints)
    return lllh


def compute_phase_affiliation(allele_matrix, haplotypes, breakpoints, prephasing, error_rate):
    """
    For each thread in each block computes the affiliation to each phase as given by the
    prephasing. Result is 3D-list, with dimensions being block id, thread inside block and phase
    of the prephasing, respectively.
    """
    ploidy = len(haplotypes)
    genpos = allele_matrix.getPositions()
    genpos_to_happos = {pos: i for i, pos in enumerate(genpos)}
    num_blocks = len(breakpoints) + 1
    block_starts = [b.position for b in breakpoints]  # list(sorted(breakpoints.keys()))
    assert block_starts == sorted(block_starts)

    # aff[b][t][p] = affinity of t-th thread in block b to p-th prephase
    aff = [[[0 for _ in range(ploidy)] for _ in range(ploidy)] for _ in range(num_blocks)]
    olp = [[[0 for _ in range(ploidy)] for _ in range(ploidy)] for _ in range(num_blocks)]
    err = [[[0 for _ in range(ploidy)] for _ in range(ploidy)] for _ in range(num_blocks)]

    # count overlaps and errors between block-threads and prephasings
    for i, pos in enumerate(prephasing.getPositions()):
        if pos not in genpos_to_happos:
            continue
        hap_pos = genpos_to_happos[pos]
        block_id = bisect_left(block_starts, hap_pos)
        for thread_id in range(ploidy):
            h_allele = haplotypes[thread_id][hap_pos]
            if h_allele < 0:
                continue
            for phase_id in range(ploidy):
                p_allele = prephasing.getAllele(phase_id, i)
                if p_allele < 0:
                    continue
                olp[block_id][thread_id][phase_id] += 1
                err[block_id][thread_id][phase_id] += 1 if h_allele != p_allele else 0

    for b in range(num_blocks):
        for t in range(ploidy):
            for p in range(ploidy):
                logprob = log(1 - error_rate) * (olp[b][t][p] - err[b][t][p])
                logprob += log(error_rate) * err[b][t][p]
                aff[b][t][p] = logprob
    return aff


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


def get_optimal_permutations(breakpoints, lllh, ploidy, affiliations):
    """
    Computes optimal permutations of haplotypes within blocks determined by breakpoints. Result
    is a list with one permutation (list) per block.
    """

    if not breakpoints:
        return [list(range(ploidy))]

    P = list(range(ploidy))
    B = list(range(len(breakpoints)))
    BE = list(range(len(breakpoints) + 1))

    # setup ILP model
    model = LpProblem(f"PermuteBlocks_p{ploidy}_b{len(breakpoints)}", LpMaximize)

    # x[b][t][h] = 1 if thread t is placed on haplotype h for block b
    x = [[[LpVariable(f"x_{b}_{t}_{h}", 0, 1, LpInteger) for h in P] for t in P] for b in BE]

    # y[b][t1][t2] = 1 if thread t1 is linked to thread t2 over breakpoint b
    y = [[[LpVariable(f"y_{b}_{t1}_{t2}", 0, 1, LpInteger) for t2 in P] for t1 in P] for b in B]

    # z[b][i] = 1 if i-th combination is used to connect haplotypes over breakpoint b
    z = [[LpVariable(f"z_{b}_{i}", 0, 1, LpInteger) for i in range(len(lllh[b]))] for b in B]

    aff_scores = []
    if affiliations is None:
        # just fix first permutation if not prephasing exists
        for t in P:
            model += x[0][t][t] == 1
    else:
        # collect all assignment-affiliation-score pairs (or force 0 if assignment not allowed)
        for b in BE:
            for t in P:
                for h in P:
                    aff_scores.append(x[b][t][h] * affiliations[b][h][t])

    # for each block there must be one-to-one-relationship between threads and haplotypes
    for i in BE:
        for j in P:
            model += sum([x[i][j][k] for k in P]) == 1
            model += sum([x[i][k][j] for k in P]) == 1

    # iff t1 and t2 both cover haplotype h on blocks b and b+1, set relationship
    for b, affected in enumerate(
        [b.haplotypes for b in breakpoints]
    ):  # enumerate(breakpoints.values()):
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
    for b, left in enumerate([b.haplotypes for b in breakpoints]):
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
    model += sum([var * weight for (var, weight) in z_weights.items()]) + sum(aff_scores)

    # solve model
    solver = PULP_CBC_CMD(msg=0)
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

    return assignments


def permute_blocks(threads, haplotypes, breakpoints, lllh, perms):

    # shuffle threads and haplotypes according to optimal assignments
    ploidy = len(haplotypes)
    threads_copy = deepcopy(threads)
    haplotypes_copy = deepcopy(haplotypes)
    ext_bp = [0] + [b.position for b in breakpoints] + [len(threads)]
    for i, (s, e) in enumerate(zip(ext_bp[:-1], ext_bp[1:])):
        for p in range(s, e):
            for t in list(range(ploidy)):
                threads[p][t] = threads_copy[p][perms[i][t]]
                haplotypes[t][p] = haplotypes_copy[perms[i][t]][p]

    # collect reorder events
    for i, bp in enumerate(breakpoints):
        # compute confidence of decision
        affected = bp.haplotypes
        assert len(lllh[i].values()) >= 2
        best = max(lllh[i].values())
        reduced = [j for j in perms[i + 1] if j in affected]
        link = tuple(affected[reduced.index(j)] for j in perms[i] if j in affected)
        bp.confidence = exp(lllh[i][link] - best) / sum([exp(v - best) for v in lllh[i].values()])
