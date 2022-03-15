import itertools as it
import logging
from collections import defaultdict
from math import log, exp
from enum import Enum

from whatshap.polyphaseutil import get_position_map

logger = logging.getLogger(__name__)


class ReorderType(Enum):
    SNP = 1
    EXIT = 2
    MULTI = 3


class ReorderEvent:
    def __init(self, position, event_type, llh_diff, haplotypes):
        self.position = position
        self.event_type = event_type
        self.confidence = 1 - exp(-llh_diff)
        self.haplotypes = haplotypes[:]


def run_reordering(
    readset,
    clustering,
    path,
    haplotypes,
    block_cut_sensitivity,
    error_rate=0.05,
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
        path, haplotypes, cwise_snps, clustering, readset, index, error_rate, window_size=20
    )

    # we can look at the sequences again to use the most likely continuation for ambiguous haplotype switches
    improve_path_on_ambiguous_switches(
        path, haplotypes, clustering, readset, index, error_rate, window_size=20
    )

    cut_positions, haploid_cuts = compute_cut_positions(
        path, block_cut_sensitivity, len(clustering)
    )

    logger.debug("Cut positions: {}".format(cut_positions))
    ploidy = len(haplotypes)
    for i in range(ploidy):
        logger.debug("Cut positions on phase {}: {}".format(i, haploid_cuts[i]))

    return (cut_positions, haploid_cuts, path, haplotypes)


def find_cluster_snps(path, haplotypes):
    """
    For each cluster, finds the positions, where it has a multiplicity of at least 2 with at least 2 different alleles.
    These positions have to be phased within the clusters.
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

    # events = []

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
                            * (error_rate ** errors)
                        )
                    score += log(likelihood)
                configs.append((perm, score))

            configs.sort(key=lambda x: -x[1])
            best_perm = configs[0][0]

            # print("   Best permutation is {} with score {}".format(best_perm, configs[0][1]))
            # print("   Second best permutation is {} with score {}".format(configs[1][0] if best_perm == configs[0][0] else configs[0][0], configs[1][1]))
            # print("   {}".format(configs))

            # switch alleles at current position, needed as input for next position(s)
            alleles = [haplotypes[best_perm[j]][pos] for j in range(len(best_perm))]
            for j in range(len(c_slots)):
                haplotypes[c_slots[j]][pos] = alleles[j]

            # events.append(ReorderEvent(pos, ReorderType.SNP, configs[0][1] - configs[1][1], set([cid])))


def improve_path_on_ambiguous_switches(
    path, haplotypes, clustering, readset, index, error_rate, window_size
):
    """
    Post processing step after the threading. If a haplotype leaves a cluster on a collapsed region, we could use
    read information to find the most likely continuation of the switching haplotypes.
    """
    if len(path) == 0:
        return []

    ploidy = len(path[0])

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
        configs = solve_single_ambiguous_site(
            haplotypes,
            inverse_perm,
            clustering,
            readset,
            index,
            i,
            h_group,
            c_group,
            error_rate,
            window_size,
        )
        best_perm = configs[0][0]

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


def solve_single_ambiguous_site(
    haplotypes,
    hap_perm,
    clustering,
    readset,
    index,
    pos,
    h_group,
    c_group,
    error_rate,
    window_size,
):
    # for every group of haplotypes coming from a collapsed cluster:
    # find the position, where each of the haplotypes joined the collapsing one
    # present_c = [corrected_path[pos - 1][j] for j in h_group]
    # next_c = [corrected_path[pos][j] for j in h_group]
    # all_c = set(present_c + next_c)

    # print("Collapsed cut on position {} on haplotypes {}:".format(pos, h_group))
    # print("   h_group = {}, present_c = {}, next_c = {}".format(h_group, present_c, next_c))

    het_pos_before = []
    het_pos_after = []

    # find the last window_size positions (starting from pos - 1),
    # which are heterozyguous among the haplotype group
    j, w = pos - 1, 0
    while w < window_size and j >= 0:
        if len(set([haplotypes[h][j] for h in h_group])) > 1:
            het_pos_before.append(j)
            w += 1
        j -= 1

    # find the next window_size positions (starting from pos),
    # which are heterozyguous among the haplotype group
    j, w = pos, 0
    het_pos_before = het_pos_before[::-1]
    while w < window_size and j < len(haplotypes[0]):
        if len(set([haplotypes[hap_perm[h]][j] for h in h_group])) > 1:
            het_pos_after.append(j)
            w += 1
        j += 1

    het_pos = het_pos_before + het_pos_after
    # print("   het_before = {}, het_after = {}".format(het_pos_before, het_pos_after))

    # store the allele of every cluster-read for the relevant positions
    rmat = dict()

    # for each position store the supporting reads
    readlist = []

    # go over all reads and determine which snp positions they are defined on
    if len(het_pos_before) > 0 and len(het_pos_after) > 0:
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
    # print("   Found {} relevant reads".format(len(readlist)))
    configs = []
    # enumerate all permutations of haplotypes at current positions
    for perm in it.permutations(range(len(h_group))):
        score = 0.0
        # for each read: determine number of errors for best fit
        for rid in readlist:
            likelihood = 0.0
            a_priori = 1 / len(h_group)
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
                likelihood += (
                    a_priori * ((1 - error_rate) ** (overlap - errors)) * (error_rate ** errors)
                )
            score += log(likelihood)
        configs.append((perm, score))

    configs.sort(key=lambda x: -x[1])
    # print("   Configs = {}".format(configs))
    return configs


def tie_break_optimal_permutations(path, pos, h_group, opt_perms):
    best_perm = opt_perms[0]
    best_score = len(path) * len(h_group) + 1
    for perm in opt_perms:
        score = 0
        for j in range(len(h_group)):
            old = h_group[perm[j]]
            new = h_group[j]
            i = pos - 1
            while i >= 0 and path[i][old] != path[pos][new]:
                i -= 1
            score += pos - 1 - i
        if score < best_score:
            best_score = score
            best_perm = perm

    return best_perm


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
