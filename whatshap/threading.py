import itertools as it
import logging
from collections import defaultdict
from .core import HaploThreader
from math import ceil, log

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
    allele_depths, normalized_depths, consensus_lists = get_allele_depths(
        readset, clustering, cov_map, ploidy, genotypes=genotypes
    )

    # compute threading through the clusters
    affine_switch_cost = ceil(compute_readlength_snp_distance_ratio(readset) / 4.0)
    path = compute_threading_path(
        readset,
        num_vars,
        cov_map,
        allele_depths,
        consensus_lists,
        ploidy,
        genotypes,
        switch_cost=4 * affine_switch_cost,
        affine_switch_cost=affine_switch_cost,
    )

    # determine haplotypes/alleles for each position
    haplotypes = compute_haplotypes(path, cov_map, consensus_lists, ploidy)

    # determine snp positions inside clusters
    cwise_snps = find_cluster_snps(path, haplotypes)

    # phase cluster snps
    haplotypes = phase_cluster_snps(path, haplotypes, cwise_snps, clustering, readset, index, window_size=5)

    # we can look at the sequences again to use the most likely continuation for ambiguous haplotype switches
    path, haplotypes = improve_path_on_ambiguous_switches(path, haplotypes, clustering, readset, index, error_rate=0.05, window_size=20)

    cut_positions, haploid_cuts = compute_cut_positions(path, block_cut_sensitivity, len(clustering))

    logger.debug("Cut positions: {}".format(cut_positions))
    for i in range(ploidy):
        logger.debug("Cut positions on phase {}: {}".format(i, haploid_cuts[i]))

    for i in range(ploidy):
        haplotypes[i] = "".join(map(lambda x: str(x) if x >= 0 else "n", haplotypes[i]))

    return (cut_positions, haploid_cuts, path, haplotypes)


def compute_readlength_snp_distance_ratio(readset):
    length = 0
    for read in readset:
        length += len(read)
    return length / len(readset)


def compute_threading_path(
    readset,
    num_vars,
    cov_map,
    allele_depths,
    consensus_lists,
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

    # run threader
    row_limit = 16 * 2 ** ploidy if ploidy > 6 else 0
    threader = HaploThreader(ploidy, switch_cost, affine_switch_cost, True, row_limit)

    path = threader.computePathsBlockwise([0], cov_map, allele_depths, consensus_lists, genotypes)
    assert len(path) == num_vars

    return path


def compute_haplotypes(path, cov_map, consensus_lists, ploidy):
    """
    Fills each haplotypes using the computed clusters and their consensus lists
    """
    num_vars = len(cov_map)
    haplotypes = [[] for _ in range(ploidy)]
    rev_cov_map = [dict() for _ in range(num_vars)]
    for i, m in enumerate(cov_map):
        for j in range(len(m)):
            rev_cov_map[i][cov_map[i][j]] = j

    for pos in range(len(path)):
        cnts = defaultdict(int)
        for i in range(ploidy):
            cid = path[pos][i]
            allele = consensus_lists[pos][rev_cov_map[pos][cid]][cnts[cid]]
            cnts[cid] += 1
            haplotypes[i].append(allele)

    return haplotypes


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


def phase_cluster_snps(path, haplotypes, cwise_snps, clustering, readset, index, window_size=10):
    """
    For clusters with multiplicity >= 2 on multiple positions: Bring the alleles of the
    precomputed haplotypes in order by using read information.
    The window size determines how many previous positions are taken into account when resolving
    the next one.
    """
    
    # sort clusters by position of first haplotype
    collapsed = sorted([(cid, sorted(snps)) for cid, snps in cwise_snps.items()], key=lambda x: x[1][0])

    # iterate cluster-wise
    for cid, snps in collapsed:
        # determine the last window_size many heterozyguous positions on the threads of the first snp
        c_slots = [j for j in range(len(path[snps[0]])) if path[snps[0]][j] == cid]
        het_pos = []
        i, w = snps[0] - 1, 0
        while w < window_size and i >= 0:
            if any([haplotypes[c_slots[0]][i] != haplotypes[c_slots[j]][i] for j in range(1, len(c_slots))]):
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
            
            best_perm = c_slots
            best_score = float("inf")
            # enumerate all permutations of haplotypes at current positions
            for perm in it.permutations(c_slots):
                score = 0
                # for each read: determine number of errors for best fit
                for rid in readlist[het_idx]:
                    if len(rmat[rid]) < 2:
                        continue
                    least_errors = window_size + 1
                    for slot in range(len(c_slots)):
                        # skip if read has different allele at current pos than tested config
                        if haplotypes[perm[slot]][pos] != rmat[rid][het_idx]:
                            continue
                        # else, compare to window_size many previous positions
                        errors = 0
                        for j in range(max(0, het_idx - window_size), het_idx):
                            if j in rmat[rid] and haplotypes[c_slots[slot]][het_pos[j]] != rmat[rid][j]:
                                errors += 1
                        least_errors = min(least_errors, errors)
                    score += least_errors
                if score < best_score:
                    best_score = score
                    best_perm = perm
            
            # switch alleles at current position, needed as input for next position(s)
            alleles = [haplotypes[best_perm[j]][pos] for j in range(len(best_perm))]
            for j in range(len(c_slots)):
                haplotypes[c_slots[j]][pos] = alleles[j]

    return haplotypes


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


def improve_path_on_ambiguous_switches(path, haplotypes, clustering, readset, index, error_rate=0.07, window_size=10):
    """
    Post processing step after the threading. If a haplotype leaves a cluster on a collapsed region, we could use
    read information to find the most likely continuation of the switching haplotypes.
    """
    if len(path) == 0:
        return []

    ploidy = len(path[0])
    corrected_path = []
    corrected_path.append(path[0])
    corrected_haplotypes = [[haplotypes[i][0]] for i in range(ploidy)]
    
    # keeps track of how the input haplotype is mapped to corrected one at current position i
    current_perm = tuple(range(ploidy))
    # the invers permutation from above, maps corrected haplotype to original one
    inverse_perm = [i for i in range(ploidy)]

    copynrs = []
    for i in range(0, len(path)):
        copynr = defaultdict(int)
        for j in range(0, ploidy):
            copynr[path[i][j]] += 1
        copynrs.append(copynr)

    for i in range(1, len(path)):
        # append next position using current permutation
        corrected_path.append([path[i][j] for j in inverse_perm])
        for j in range(len(inverse_perm)):
            corrected_haplotypes[j].append(haplotypes[inverse_perm[j]][i])
        
        # iterate over present cluster ids
        changed = set()
        for c_id in copynrs[i]:
            if copynrs[i - 1][c_id] >= 2:
                # for all collapsed clusters: find haplotypes, which go through and check whether one of them exits
                outgoing_c = False
                affected = []
                for j in range(ploidy):
                    if corrected_path[i - 1][j] == c_id:
                        affected.append(j)
                        if corrected_path[i][j] != c_id:
                            outgoing_c = True
                # if haplotypes leaves collapsed cluster, all other might be equally suited, so add them
                if outgoing_c:
                    for h in affected:
                        changed.add(h)

        for j in range(ploidy):
            if corrected_path[i - 1][j] != corrected_path[i][j]:
                changed.add(j)

        # if only haplotype switched cluster and only from non-collapsed cluster: skip
        if len(changed) < 2:
            continue

        # compute the optimal permutation of changing haplotypes, according to read information
        h_group = list(changed)
        configs = solve_single_ambiguous_site(corrected_path, haplotypes, corrected_haplotypes, inverse_perm, clustering, readset, index, i, h_group, error_rate, window_size)
        
        # if multiple optimal configs exist: let the most recently haplotypes switch
        opt_perms = [config[0] for config in configs if config[1] == configs[0][1]]
        if len(opt_perms) > 1:
            best_perm = tie_break_optimal_permutations(corrected_path, i, h_group, opt_perms)
        else:
            best_perm = configs[0][0]

        #print("   Best permutation is {} with score {}".format(best_perm, configs[0][1]))
        #print("   Second best permutation is {} with score {}".format(configs[1][0] if best_perm == configs[0][0] else configs[0][0], configs[1][1]))
            
        # apply local best permutation to current global permutation
        current_perm_copy = list(current_perm)
        for j in range(len(h_group)):
            current_perm_copy[inverse_perm[h_group[j]]] = h_group[best_perm[j]]
        current_perm = tuple(current_perm_copy)
        for j in range(ploidy):
            inverse_perm[current_perm[j]] = j

        # correct current position according to updated permutation
        for j in range(len(inverse_perm)):
            corrected_path[i][j] = path[i][inverse_perm[j]]
            corrected_haplotypes[j][i] = haplotypes[inverse_perm[j]][i]

    assert len(path) == len(corrected_path)
    assert len(haplotypes[0]) == len(corrected_haplotypes[0])
    return corrected_path, corrected_haplotypes


def solve_single_ambiguous_site(corrected_path, haplotypes, corrected_haplotypes, hap_perm, clustering, readset, index, pos, h_group, error_rate, window_size):
    # for every group of haplotypes coming from a collapsed cluster:
    # find the position, where each of the cluster joined the collapsing one (introduced in new07)
    present_c = [corrected_path[pos - 1][j] for j in h_group]
    next_c = [corrected_path[pos][j] for j in h_group]
    all_c = set(present_c + next_c)

    #print("Collapsed cut on position {} on haplotypes {}:".format(pos, h_group))
    #print("   h_group = {}, present_c = {}, next_c = {}".format(h_group, present_c, next_c))

    het_pos_before = []
    het_pos_after = []

    j, w = pos - 1, 0
    while w < window_size and j >= 0:
        if any([corrected_haplotypes[h_group[0]][j] != corrected_haplotypes[h_group[k]][j] for k in range(1, len(h_group))]):
            het_pos_before.append(j)
            w += 1
        j -= 1
    j, w = pos, 0
    het_pos_before = het_pos_before[::-1]
    while w < window_size and j < len(haplotypes[0]):
        if any([haplotypes[hap_perm[h_group[0]]][j] != haplotypes[hap_perm[h_group[k]]][j] for k in range(1, len(h_group))]):
            het_pos_after.append(j)
            w += 1
        j += 1

    het_pos = het_pos_before + het_pos_after
    #print("   het_before = {}, het_after = {}".format(het_pos_before, het_pos_after))

    # store the allele of every cluster-read for the relevant positions
    rmat = dict()

    # for each position store the supporting reads
    readlist = []

    # go over all reads and determine which snp positions they are defined on
    if len(het_pos_before) > 0 and len(het_pos_after) > 0:
        for cid in all_c:
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
    #print("   Found {} relevant reads".format(len(readlist)))
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
                for j in het_pos_before:
                    if j in rmat[rid] and corrected_haplotypes[h_group[perm[slot]]][j] != rmat[rid][j]:
                        errors += 1
                for j in het_pos_after:
                    if j in rmat[rid] and haplotypes[hap_perm[h_group[slot]]][j] != rmat[rid][j]:
                        errors += 1
                likelihood += a_priori * ((1 - error_rate)**(overlap - errors)) * (error_rate**errors)
            score += log(likelihood)
        configs.append((perm, score))
            
    configs.sort(key=lambda x: -x[1])
    #print("   Configs = {}".format(configs))
    return configs


def tie_break_optimal_permutations(corrected_path, pos, h_group, opt_perms):
    best_perm = opt_perms[0]
    best_score = len(corrected_path) * len(h_group) + 1
    for perm in opt_perms:
        score = 0
        for j in range(len(h_group)):
            old = h_group[perm[j]]
            new = h_group[j]
            i = pos - 1
            while i >= 0 and corrected_path[i][old] != corrected_path[pos][new]:
                i -= 1
            score += pos - 1 - i
        if score < best_score:
            best_score = score
            best_perm = perm
            
    return best_perm


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

    # re-add clusters with 1-position-deletion
    for pos in range(1, len(cov_map) - 1):
        for cid in cov_map[pos - 1]:
            if cid in cov_map[pos + 1] and cid not in cov_map[pos]:
                if len(cov_map[pos]) > cut_off:
                    break
                cov_map[pos].append(cid)

    return cov_map


def get_coverage(readset, clustering, pos_index):
    """
    Returns a list, which for every position contains a dictionary, mapping a cluster id to
    a relative coverage on this position.
    """
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


def get_coverage_absolute(readset, clustering, pos_index):
    """
    Returns a list, which for every position contains a dictionary, mapping a cluster id to
    an absolute coverage on this position.
    """
    num_vars = len(pos_index)
    num_clusters = len(clustering)
    coverage = [defaultdict(int) for pos in range(num_vars)]
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


def get_allele_depths(readset, clustering, cov_map, ploidy, genotypes=None):
    """
    Returns a list, which for every position contains a list (representing the clusters) of dictionaries containing the allele depths.
    Additionally computes a consensus list per position per cluster, such that the first k elements represent the alleles of this
    cluster, if it was selected with multiplicity k.
    ad[pos][c_id][al] = number of reads in cluster cov_map[pos][c_id] having allele al at position pos
    Indices are local, i.e. the i-th entry of ad is the entry for the i-th position that occurs in the readset.

    If genotypes are provided, the allele depths are normalized per position, such that the ratio between alleles is (roughly) the
    same as according to the genotypes.
    """
    # Map genome positions to [0,l)
    index = {}
    rev_index = []
    num_vars = 0
    for position in readset.get_positions():
        index[position] = num_vars
        rev_index.append(position)
        num_vars += 1

    # stores allele depth per position and cluster
    ad = [[dict() for c_id in cov_map[pos]] for pos in range(num_vars)]
    nad = [[dict() for c_id in cov_map[pos]] for pos in range(num_vars)]
    # store depths per position (over all clusters)
    ad_per_pos = [dict() for pos in range(num_vars)]

    # Create reverse map of the used clusters for every position
    rev_cov_map = [dict() for _ in range(num_vars)]
    for i, m in enumerate(cov_map):
        for j in range(len(m)):
            rev_cov_map[i][cov_map[i][j]] = j

    # count alleles
    for c_id, cluster in enumerate(clustering):
        for read in cluster:
            for var in readset[read]:
                pos = index[var.position]
                if c_id in rev_cov_map[pos]:
                    al = var.allele
                    if al not in ad[pos][rev_cov_map[pos][c_id]]:
                        ad[pos][rev_cov_map[pos][c_id]][al] = 0
                    ad[pos][rev_cov_map[pos][c_id]][al] += 1
                    if al not in ad_per_pos[pos]:
                        ad_per_pos[pos][al] = 0
                    ad_per_pos[pos][al] += 1

    # normalize if genotypes provided
    if genotypes:
        assert len(genotypes) == len(ad)
        for pos in range(len(ad)):
            coverage = sum(cnt for cnt in ad_per_pos[pos].values())
            exp_ad = {al: genotypes[pos][al] * coverage / ploidy for al in ad_per_pos[pos]}
            for cid in range(len(ad[pos])):
                for al in ad[pos][cid]:
                    nad[pos][cid][al] = round(ad[pos][cid][al] * exp_ad[al] / ad_per_pos[pos][al])

    # compute allele lists
    cons_lists = [[[] for c_id in cov_map[pos]] for pos in range(num_vars)]
    for pos in range(num_vars):
        for c_id in range(len(cov_map[pos])):
            cnts = defaultdict(int)
            for i in range(ploidy):
                max_cnt = 0
                max_al = 0
                for al in nad[pos][c_id]:
                    cnt = nad[pos][c_id][al] / (1 + cnts[al])
                    if cnt > max_cnt:
                        max_cnt = cnt
                        max_al = al
                cons_lists[pos][c_id].append(max_al)
                cnts[max_al] += 1

    return ad, nad, cons_lists


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
        newdict = defaultdict(int)
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
