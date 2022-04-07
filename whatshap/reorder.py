import itertools as it
import logging
from collections import defaultdict
from math import log, exp
from enum import Enum

from whatshap.polyphaseutil import get_position_map

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
    readset,
    clustering,
    path,
    haplotypes,
    block_cut_sensitivity,
    error_rate=0.025,
):
    """
    Main method for the reordering stage of the polyploid phasing algorithm. Input:

    readset -- The fragment matrix to phase
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

    # compute auxiliary data
    index, rev_index = get_position_map(readset)

    # determine snp positions inside clusters
    cwise_snps = find_cluster_snps(path, haplotypes)

    # phase cluster snps
    phase_cluster_snps(
        path, haplotypes, cwise_snps, clustering, readset, index, error_rate, window_size=32
    )

    # select most likely continuation for ambiguous haplotype switches
    events = resolve_ambiguous_switches(
        path, haplotypes, clustering, readset, index, error_rate, window_size=32
    )

    cut_positions, haploid_cuts = compute_cut_positions(path, block_cut_sensitivity, events)

    logger.debug("Cut positions: {}".format(cut_positions))
    ploidy = len(haplotypes)
    for i in range(ploidy):
        logger.debug("Cut positions on phase {}: {}".format(i, haploid_cuts[i]))

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
    path, haplotypes, cwise_snps, clustering, readset, index, error_rate, window_size
):
    """
    For clusters with multiplicity >= 2 on multiple positions: Bring the alleles of the
    precomputed haplotypes in order by using read information.
    The window size determines how many previous positions are taken into account when resolving
    the next one.
    """
    # sort clusters by position of first haplotype
    collapsed = sorted(
        [(cid, sorted(snps)) for cid, snps in cwise_snps.items()], key=lambda x: x[1][0]
    )

    # iterate cluster-wise
    for cid, snps in collapsed:
        # determine the last window_size many heterozyguous positions on the threads of the first snp
        c_slots = [j for j in range(len(path[snps[0]])) if path[snps[0]][j] == cid]
        het_pos = []
        i, w = snps[0] - 1, 0
        while w < window_size and i >= 0:
            if len(set([haplotypes[s][i] for s in c_slots])) > 1:
                het_pos.append(i)
                w += 1
            i -= 1
        het_pos = het_pos[::-1] + snps
        # w = number of preceding het positions, can be lower than window_size close to chr start

        # store the allele of every cluster-read for the relevant positions
        rmat = dict()

        # for each position store the supporting reads
        readlist = defaultdict(list)

        # go over all reads and determine which snp positions they are defined on
        for rid in clustering[cid]:
            read = readset[rid]
            rmat[rid] = dict()
            i, j = 0, 0
            while i < len(read) and j < len(het_pos):
                read_pos = index[read[i].position]
                if read_pos == het_pos[j]:
                    rmat[rid][j] = read[i].allele
                    readlist[j].append(rid)
                    i += 1
                    j += 1
                elif read_pos < het_pos[j]:
                    i += 1
                else:
                    j += 1

        # reconstruct cluster phasing
        for i in range(len(snps)):
            pos = snps[i]
            het_idx = i + w
            c_slots = [j for j in range(len(path[pos])) if path[pos][j] == cid]

            configs = []
            # enumerate all permutations of haplotypes at current positions
            for perm in it.permutations(c_slots):
                score = 0
                # for each read: determine number of errors for best fit
                for rid in readlist[het_idx]:
                    if len(rmat[rid]) < 2:
                        continue
                    likelihood = 0.0
                    a_priori = 1 / len(c_slots)
                    for slot in range(len(c_slots)):
                        overlap, errors = 0, 0
                        # skip if read has different allele at current pos than tested config
                        if haplotypes[perm[slot]][pos] != rmat[rid][het_idx]:
                            errors += 1
                        # else, compare to window_size many previous positions
                        for j in range(max(0, het_idx - window_size), het_idx):
                            if j in rmat[rid]:
                                overlap += 1
                                if haplotypes[c_slots[slot]][het_pos[j]] != rmat[rid][j]:
                                    errors += 1
                        likelihood += (
                            a_priori
                            * ((1 - error_rate) ** (overlap - errors))
                            * (error_rate**errors)
                        )
                    score += log(likelihood)
                configs.append((perm, score))

            configs.sort(key=lambda x: -x[1])
            best_perm = configs[0][0]

            # switch alleles at current position, needed as input for next position(s)
            alleles = [haplotypes[best_perm[j]][pos] for j in range(len(best_perm))]
            for j in range(len(c_slots)):
                haplotypes[c_slots[j]][pos] = alleles[j]


def resolve_ambiguous_switches(
    path, haplotypes, clustering, readset, index, error_rate, window_size
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
        changed_h = set([j for j in range(ploidy) if path[i - 1][j] != path[i][inverse_perm[j]]])
        c_group = set([path[i - 1][j] for j in changed_h])
        h_group = sorted([j for j in range(ploidy) if path[i - 1][j] in c_group])

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
                readset,
                index,
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
    readset,
    index,
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
        while len(set([path[j][h] for h in h_group])) < 2 and j > 0:
            j -= 1

    while len(het_pos_before) < window_size and j >= 0:
        if len(set([haplotypes[h][j] for h in h_group])) > 1:
            het_pos_before.append(j)
        j -= 1
    het_pos_before = het_pos_before[::-1]
    # same for the next window_size positions (starting from pos)
    j = pos
    while len(het_pos_after) < window_size and j < len(haplotypes[0]):
        if len(set([haplotypes[hap_perm[h]][j] for h in h_group])) > 1:
            het_pos_after.append(j)
        j += 1
    het_pos = het_pos_before + het_pos_after

    # store the allele of every cluster-read for the relevant positions
    rmat = dict()

    # for each position store the supporting reads
    readlist = []

    # go over all reads and determine which snp positions they are defined on
    if len(het_pos_before) >= 2 and len(het_pos_after) >= 2:
        for cid in c_group:
            for rid in clustering[cid]:
                read = readset[rid]
                if len(read) < 2:
                    continue
                if index[read[0].position] > het_pos_before[-1]:
                    continue
                if index[read[-1].position] < het_pos_after[0]:
                    continue
                rmat[rid] = dict()
                j, k, b, t = 0, 0, 0, 0
                while j < len(read) and k < len(het_pos):
                    read_pos = index[read[j].position]
                    if read_pos == het_pos[k]:
                        rmat[rid][het_pos[k]] = read[j].allele
                        b += 1 if k < len(het_pos_before) else 0
                        t, j, k = t + 1, j + 1, k + 1
                    elif read_pos < het_pos[k]:
                        j += 1
                    else:
                        k += 1

                if b > 0 and t > b:
                    readlist.append(rid)

    # reconstruct cluster phasing
    configs = []
    # enumerate all permutations of haplotypes at current positions
    for perm in it.permutations(range(len(h_group))):
        score = 0.0
        num_factors = 0
        # for each read: determine number of errors for best fit
        for rid in readlist:
            likelihood = 0.0
            for slot in range(len(h_group)):
                overlap = len(rmat[rid])
                errors = 0
                for p in rmat[rid]:
                    if p < het_pos_after[0]:
                        cmp = haplotypes[h_group[perm[slot]]][p]
                    else:
                        cmp = haplotypes[hap_perm[h_group[slot]]][p]
                    if cmp != rmat[rid][p]:
                        errors += 1
                likelihood = max(
                    likelihood, ((1 - error_rate) ** (overlap - errors)) * (error_rate**errors)
                )
            num_factors += overlap
            score += log(likelihood)
        configs.append((perm, score / num_factors if num_factors > 0 else 0))

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
