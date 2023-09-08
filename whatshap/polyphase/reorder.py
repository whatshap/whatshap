import itertools as it
import logging
from collections import defaultdict
from bisect import bisect_right
from typing import List
from math import log, exp
from functools import reduce
from operator import mul
from copy import deepcopy
from pulp import LpProblem, LpVariable, LpMaximize, LpInteger

from whatshap.polyphase import PolyphaseBlockResult, PhaseBreakpoint, get_ilp_solver

logger = logging.getLogger(__name__)


def find_subinstances(allele_matrix, clustering, threads, haplotypes):
    """
    Scans clusters for heterozygous positions, i.e. the cluster contains at least 2 threads,
    but the alleles are different. Returns a list of triplets, containing of a cluster id,
    a set of affected threads and a submatrix with only the heterozygous positions and reads
    covering them.
    Each triplet represents a sub-instance to phase to determine allele order inside the cluster.
    A cluster can be present in multiple sub-instances if the set of threads inside the cluster
    changed. In this case, the results of both regions cannot be combined in a meaningful way.
    """

    # iterate over positions, keep a ploidy and list of found positions for each cluster
    cwise_snps = defaultdict(list)
    last_thread_set = defaultdict(list)
    collapsed = []
    for pos in range(len(threads)):
        clusters = set()
        alleles = defaultdict(set)
        thread_set = defaultdict(list)
        for i, cid in enumerate(threads[pos]):
            clusters.add(cid)
            alleles[cid].add(haplotypes[i][pos])
            thread_set[cid].append(i)
        for cid in clusters:
            if len(alleles[cid]) >= 2:
                # if thread-set changed: write old list to results and start new one
                if last_thread_set[cid] != thread_set[cid]:
                    if cwise_snps[cid]:
                        collapsed.append((cid, last_thread_set[cid], cwise_snps[cid]))
                    last_thread_set[cid] = thread_set[cid]
                    cwise_snps[cid] = []
                cwise_snps[cid].append(pos)

    # write remaining lists into collapsed
    for cid, snps in cwise_snps.items():
        if snps:
            assert len(last_thread_set[cid]) > 0
            collapsed.append((cid, last_thread_set[cid], snps))

    sub_instances = []
    num_vars = len(allele_matrix.getPositions())
    ploidy = len(haplotypes)
    for cid, thread_set, snps in collapsed:
        if len(snps) == num_vars and len(thread_set) == ploidy:
            continue
        subm = allele_matrix.extractSubMatrix(snps, clustering[cid], True)
        assert len(subm.getPositions()) > 0
        if len(subm) > 0:
            sub_instances.append((cid, thread_set, subm))

    return sub_instances


# def integrate_sub_results(allele_matrix, sub_instances, sub_results, threads, haplotypes):
def integrate_sub_results(
    allele_matrix,
    sub_instances,
    sub_results: List[PolyphaseBlockResult],
    threads: List[List[int]],
    haplotypes,
):
    """
    Does two things:
    1. Update haplotype strings inside the collapsed regions according to the solved sub-instances
    2. Collect breakpoints of global threading and integrate breakpoints of sub-instances
    """
    breakpoints = find_breakpoints(threads)
    for (cid, thread_set, subm), res in zip(sub_instances, sub_results):
        snps = [allele_matrix.globalToLocal(gpos) for gpos in subm.getPositions()]
        assert all([0 <= pos < allele_matrix.getNumPositions() for pos in snps])

        # reorder haplotype alleles according to subresults, but only on SNP positions
        for i, pos in enumerate(snps):
            for j, hap in enumerate(thread_set):
                haplotypes[hap][pos] = res.haplotypes[j][i]

        # copy breakpoints, but map positions and haplotype ids from sub-instance to global one
        for bp in res.breakpoints:
            pos = allele_matrix.globalToLocal(subm.localToGlobal(bp.position))
            haps = [thread_set[i] for i in bp.haplotypes]
            breakpoints.append(PhaseBreakpoint(pos, haps, bp.confidence))

    # Join duplicate breakpoints
    breakpoints.sort(key=lambda x: x.position)
    i = 0
    while i < len(breakpoints):
        j = i + 1
        while j < len(breakpoints) and breakpoints[i].position == breakpoints[j].position:
            j += 1
        if i + 1 == j:
            i += 1
            continue
        # breakpoints[i:j] have same position with |j - i| >= 2
        haps = sorted(list({h for k in range(i, j) for h in breakpoints[k].haplotypes}))
        conf = reduce(mul, [breakpoints[k].confidence for k in range(i, j)])
        breakpoints[i].haplotypes = haps
        breakpoints[i].confidence = conf
        del breakpoints[i + 1 : j]
        assert i + 1 >= len(breakpoints) or breakpoints[i].position != breakpoints[i + 1]
        i += 1

    return breakpoints


def run_reordering(
    allele_matrix,
    clustering,
    threads,
    haplotypes,
    breakpoints,
    prephasing,
    error_rate=0.07,
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
    ploidy = len(haplotypes)
    perms = get_optimal_assignments(breakpoints, lllh, ploidy, aff)
    permute_blocks(threads, haplotypes, breakpoints, lllh, perms)


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
    dimension contains k! entries, representing the linkage likelihoods for every possible
    permutation of the k affected haplotypes.
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
        if pos > 0:
            affected_clusts = affected_clusts.union({threads[pos - 1][h] for h in affected})
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
                for j, a in read:
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
            left_h = list(affected)
            right_h = [perm[affected.index(i)] for i in affected]
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
    prephasing_pos = prephasing.getPositions()
    phaseblock_starts = sorted(list({prephasing.getFirstPos(i) for i in range(len(prephasing))}))
    phaseblock_starts.append(len(prephasing_pos))
    for phb, (start, end) in enumerate(zip(phaseblock_starts[:-1], phaseblock_starts[1:])):
        for i in range(start, end):
            pos = prephasing_pos[i]
            if pos not in genpos_to_happos:
                continue
            hap_pos = genpos_to_happos[pos]
            block_id = bisect_right(block_starts, hap_pos)
            for thread_id in range(ploidy):
                h_allele = haplotypes[thread_id][hap_pos]
                if h_allele < 0:
                    continue
                for phase_id in range(phb * ploidy, (phb + 1) * ploidy):
                    # accessing read set (of ploidy many reads) for current phase block
                    p_allele = prephasing.getAllele(phase_id, i)
                    if p_allele < 0:
                        continue
                    olp[block_id][thread_id][phase_id % ploidy] += 1
                    err[block_id][thread_id][phase_id % ploidy] += 1 if h_allele != p_allele else 0

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


def get_optimal_assignments(breakpoints, lllh, ploidy, affiliations):
    """
    Computes optimal permutations of haplotypes within blocks determined by breakpoints. Result
    is a list with one permutation (list) per block.
    """
    P = list(range(ploidy))
    B = list(range(len(breakpoints)))
    BE = list(range(len(breakpoints) + 1))

    if not breakpoints:
        return [list(range(ploidy))]

    # Only solve ILP if pre-phasing affiliations are given. Otherwise, take local optima
    if not affiliations:
        assignments = [[i for i in P] for _ in BE]
        for b in B:
            for i in P:
                assignments[b + 1][i] = assignments[b][i]
            perm = max(lllh[b], key=lllh[b].get)
            affected = sorted(perm)
            for left, right in zip(affected, perm):
                assignments[b + 1][assignments[b].index(left)] = right

        return assignments

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
    solver = solver = get_ilp_solver()
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
