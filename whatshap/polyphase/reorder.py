import itertools as it
import logging
from collections import defaultdict
from math import log, exp
from enum import Enum

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
    path,
    haplotypes,
    block_cut_sensitivity,
    error_rate=0.025,
):
    """
    Main method for the reordering stage of the polyploid phasing algorithm. Input:

    allele_matrix -- The fragment matrix to phase
    clustering -- A list of clusters. Each cluster is a list of read ids, indicating which reads it
                  contains. Every read can only be present in one cluster.
    path -- The computed cluster paths from the threading stage
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
    cwise_snps = find_cluster_snps(path, haplotypes)

    # phase cluster snps
    phase_cluster_snps(
        path, haplotypes, cwise_snps, clustering, allele_matrix, error_rate, window_size=32
    )

    # select most likely continuation for ambiguous haplotype switches
    events = resolve_ambiguous_switches(
        path, haplotypes, clustering, allele_matrix, error_rate, window_size=32
    )

    cut_positions, haploid_cuts = compute_cut_positions(path, block_cut_sensitivity, events)

    logger.debug(f"Cut positions: {cut_positions}")
    ploidy = len(haplotypes)
    for i in range(ploidy):
        logger.debug(f"Cut positions on phase {i}: {haploid_cuts[i]}")

    return (cut_positions, haploid_cuts, path, haplotypes)


def find_cluster_snps(path, haplotypes):
    """
    For each cluster, finds the positions, where it has a multiplicity of at least 2 with at least
    2 different alleles. These positions have to be phased within the clusters.
    """
    cwise_snps = defaultdict(list)
    for pos in range(len(path)):
        clusters = set()
        alleles = defaultdict(set)
        for i, cid in enumerate(path[pos]):
            clusters.add(cid)
            alleles[cid].add(haplotypes[i][pos])
        for cid in clusters:
            if len(alleles[cid]) >= 2:
                cwise_snps[cid].append(pos)

    return cwise_snps


def phase_cluster_snps(
    path, haplotypes, cwise_snps, clustering, allele_matrix, error_rate, window_size
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
        c_slots = [j for j in range(len(path[snps[0]])) if path[snps[0]][j] == cid]
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
            c_slots = [j for j in range(len(path[pos])) if path[pos][j] == cid]

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


def resolve_ambiguous_switches(
    path, haplotypes, clustering, allele_matrix, error_rate, window_size
):
    """
    Post processing step after the threading. If a haplotype leaves a cluster on a collapsed region,
    we could use read information to find the most likely continuation of the switching haplotypes.
    """
    if len(path) == 0:
        return []

    ploidy = len(path[0])
    events = []
    event_type = ReorderType.NONE

    # maps the haplotype slot of the corrected version to the original (and in reverse)
    current_perm, inverse_perm = list(range(ploidy)), list(range(ploidy))

    for i in range(1, len(path)):
        # find slots affected by cluster changes
        changed_h = {j for j in range(ploidy) if path[i - 1][j] != path[i][inverse_perm[j]]}
        c_group = {path[i - 1][j] for j in changed_h}
        h_group = sorted(j for j in range(ploidy) if path[i - 1][j] in c_group)

        # if only one slot involved -> just write current permutation to haplotype object
        if len(h_group) < 2:
            path[i] = [path[i][j] for j in inverse_perm]
            reord = [haplotypes[inverse_perm[j]][i] for j in range(ploidy)]
            for j in range(ploidy):
                haplotypes[j][i] = reord[j]
            continue

        # compute the optimal permutation of changing haplotypes, according to read information
        for b in [True, False]:
            configs = solve_single_ambiguous_site(
                haplotypes,
                inverse_perm,
                clustering,
                allele_matrix,
                i,
                path,
                h_group,
                c_group,
                b,
                error_rate,
                window_size,
            )
            if len(c_group) == 1:
                event_type = ReorderType.COLLAPSED_PRE if b else ReorderType.COLLAPSED
            else:
                event_type = ReorderType.MULTI
            best_perm = configs[0][0]
            if configs[0][1] > configs[1][1] + 0.5:
                break
        event = ReorderEvent(i, event_type, exp(configs[0][1]), exp(configs[1][1]), h_group)
        events.append(event)

        # apply local best permutation to current global permutation
        current_perm_copy = list(current_perm)
        for j in range(len(h_group)):
            current_perm_copy[inverse_perm[h_group[j]]] = h_group[best_perm[j]]
        current_perm = current_perm_copy
        for j in range(ploidy):
            inverse_perm[current_perm[j]] = j

        # correct current position according to updated permutation
        path[i] = [path[i][j] for j in inverse_perm]
        reord = [haplotypes[inverse_perm[j]][i] for j in range(ploidy)]
        for j in range(ploidy):
            haplotypes[j][i] = reord[j]

    return events


def solve_single_ambiguous_site(
    haplotypes,
    hap_perm,
    clustering,
    allele_matrix,
    pos,
    path,
    h_group,
    c_group,
    exclude_collapsed_cluster,
    error_rate,
    window_size,
):
    # find the last window_size positions (starting from pos - 1),
    # which are heterozyguous among the haplotype group
    het_pos_before, het_pos_after = [], []
    j = pos - 1
    if exclude_collapsed_cluster:
        while len({path[j][h] for h in h_group}) < 2 and j > 0:
            j -= 1

    while len(het_pos_before) < window_size and j >= 0:
        if len({haplotypes[h][j] for h in h_group}) > 1:
            het_pos_before.append(j)
        j -= 1
    het_pos_before = het_pos_before[::-1]
    # same for the next window_size positions (starting from pos)
    j = pos
    while len(het_pos_after) < window_size and j < len(haplotypes[0]):
        if len({haplotypes[hap_perm[h]][j] for h in h_group}) > 1:
            het_pos_after.append(j)
        j += 1
    het_pos = het_pos_before + het_pos_after

    # use submatrix of relevant reads and positions
    rids = filter(
        lambda r: allele_matrix.getFirstPos(r) < pos <= allele_matrix.getLastPos(r),
        [r for cid in c_group for r in clustering[cid]],
    )
    submatrix = allele_matrix.extractSubMatrix(het_pos, list(rids), True)

    # reconstruct cluster phasing
    configs = []
    # enumerate all permutations of haplotypes at current positions
    for perm in it.permutations(range(len(h_group))):
        score = 0.0
        num_factors = 0
        # for each read: determine number of errors for best fit
        for rid in range(len(submatrix)):
            likelihood = 0.0
            read = submatrix.getRead(rid)
            overlap = len(read)
            for slot in range(len(h_group)):
                errors = 0
                for (j, a) in read:
                    p = het_pos[j]
                    if j < len(het_pos_before):
                        cmp = haplotypes[h_group[perm[slot]]][p]
                    else:
                        cmp = haplotypes[hap_perm[h_group[slot]]][p]
                    if cmp != a:
                        errors += 1
                likelihood = max(
                    likelihood, ((1 - error_rate) ** (overlap - errors)) * (error_rate**errors)
                )
            num_factors += overlap
            score += log(likelihood)
        configs.append((perm, score / num_factors if num_factors > 0 else 0))

    del submatrix
    configs.sort(key=lambda x: -x[1])
    return configs


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
