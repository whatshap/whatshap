import itertools
import logging
from collections import defaultdict
from .solver import HaploThreader
from math import ceil, log
from scipy.stats import binom

logger = logging.getLogger(__name__)


def run_threading(
    allele_matrix,
    clustering,
    ploidy,
    genotypes,
    distrust_genotypes=False,
    max_cluster_gap=10,
    error_rate=0.05,
):
    """
    Main method for the threading stage of the polyploid phasing algorithm. Taken inputs:

    allele_matrix -- The fragment matrix to phase
    clustering -- A list of clusters. Each cluster is a list of read ids, indicating which reads it
                  contains. Every read can only be present in one cluster.
    ploidy -- Number of haplotypes to phase
    genotypes -- Desired genotype as dictionary (allele -> multiplicity), as list over positions.
    distrust_genotypes -- If False, forces computed genotypes into given ones with least amount of
                          changes. Otherwise, genotypes will only be used to fill missing alleles.
    max_cluster_gap -- Number of allowed consecutive zero-coverage positions for a cluster to still
                       be considered for the threading
    error_rate -- Estimation of allele error rate for the reads. Only relevant when forcing genotypes

    """

    # compute auxiliary data
    num_vars = allele_matrix.getNumPositions()
    allele_depths, cons_lists = get_allele_depths(allele_matrix, clustering, ploidy)
    cov_map = select_clusters(allele_depths, ploidy, max_cluster_gap)

    # compute threading through the clusters
    affine_switch_cost = ceil(compute_readlength_snp_distance_ratio(allele_matrix) / 1.0)
    paths = compute_threading_path(
        cov_map,
        allele_depths,
        ploidy,
        switch_cost=4 * affine_switch_cost,
        affine_switch_cost=affine_switch_cost,
        max_cluster_gap=max_cluster_gap,
    )
    assert len(paths) == num_vars

    # determine haplotypes/alleles for each position
    haplotypes = compute_haplotypes(paths, cons_lists, ploidy)

    # enforce genotypes
    if not distrust_genotypes:
        haplotypes = force_genotypes(
            paths, haplotypes, genotypes, cov_map, allele_depths, error_rate
        )

    return (paths, haplotypes)


def compute_readlength_snp_distance_ratio(allele_matrix):
    length = 0
    for read in allele_matrix:
        length += len(read)
    return length / len(allele_matrix)


def compute_threading_path(
    cov_map,
    allele_depths,
    ploidy,
    switch_cost=32.0,
    affine_switch_cost=8.0,
    max_cluster_gap=10,
):
    """
    Runs the threading algorithm for the haplotypes using the given costs for switches. The normal switch cost is the
    cost per haplotype, that switches clusters from one position to another (non matching coverage on single position
    is always cost 1.0). The affine switch cost is an additional offset for every position, where a switch occurs.
    These additional costs encourage the threading algorithm to summarize multiple switches on consecutive positions
    into one big switch.
    """

    logger.debug(f"Computing threading paths with switch cost {switch_cost} ..")

    # run threader
    row_limit = 16 * 2**ploidy if ploidy > 6 else 0
    threader = HaploThreader(ploidy, switch_cost, affine_switch_cost, max_cluster_gap, row_limit)

    path = threader.computePathsBlockwise([0], cov_map, allele_depths)

    return path


def compute_haplotypes(path, consensus_lists, ploidy):
    """
    Fills each haplotypes using the computed clusters and their consensus lists
    """
    haplotypes = [[] for _ in range(ploidy)]

    for pos in range(len(path)):
        cnts = defaultdict(int)
        for i in range(ploidy):
            cid = path[pos][i]
            if cid in consensus_lists[pos]:
                allele = consensus_lists[pos][cid][cnts[cid]]
            else:
                allele = -1
            cnts[cid] += 1
            haplotypes[i].append(allele)

    return haplotypes


def force_genotypes(path, haplotypes, genotypes, cov_map, allele_depths, error_rate):
    num_vars = len(path)
    for pos in range(num_vars):
        # count allele occurences
        alleles = {a for a in genotypes[pos]}
        present = defaultdict(int)
        for h in haplotypes:
            present[h[pos]] += 1
            alleles.add(h[pos])

        # skip positions with undetermined allele
        if -1 in present:
            continue

        # detect abundances and shortages
        abundant_alleles, lacking_alleles = dict(), dict()
        alleles_to_insert, affected_positions = [], []
        for a in alleles:
            if a not in genotypes[pos]:
                genotypes[pos][a] = 0
            diff = present[a] - genotypes[pos][a]
            if diff > 0:
                # slots with abundant allele could change, in total we want the lower count
                abundant_alleles[a] = diff
                alleles_to_insert += [a for _ in range(genotypes[pos][a])]
                for p in range(len(path[pos])):
                    if haplotypes[p][pos] == a:
                        affected_positions.append(p)
            elif diff < 0:
                # slots with short alleles can stay, we need to insert the missing amount
                lacking_alleles[a] = -diff
                alleles_to_insert += [a for _ in range(-diff)]

        affected_positions.sort()
        alleles_to_insert.sort()

        if len(abundant_alleles) == 0:
            continue

        # for all permutations, pick the one best fitting to allele depths of clusters
        clusts = cov_map[pos]
        given_config = [haplotypes[h][pos] for h in range(len(haplotypes))]
        best_config = given_config
        best_likelihood = -float("inf")
        for perm in set(list(itertools.permutations(alleles_to_insert))):
            # build next config (given + affected slots permuted)
            newconfig = given_config[:]
            for i in range(len(perm)):
                newconfig[affected_positions[i]] = perm[i]

            # compute likelihood how well newconfig explains the observed allele depths of the clusters
            log_likelihood = 0

            # compute allele depth fraction for each cluster
            for clust in clusts:
                allele_mult = {a: 0 for a in alleles}
                clust_mult = 0
                for slot in range(len(path[pos])):
                    if path[pos][slot] == clust:
                        allele_mult[newconfig[slot]] += 1
                        clust_mult += 1
                if clust_mult > 0:
                    # compute likelihood that observed allele depth was produced by assigned alleles per cluster
                    total_depth = sum(allele_depths[pos][clust].values())
                    for a in alleles:
                        allele_mult[a] /= clust_mult
                        allele_mult[a] = (
                            allele_mult[a] * (1 - error_rate) + (1 - allele_mult[a]) * error_rate
                        )
                        observed_depth = (
                            0
                            if a not in allele_depths[pos][clust]
                            else allele_depths[pos][clust][a]
                        )
                        prob = binom.pmf(observed_depth, total_depth, allele_mult[a])
                        log_likelihood += log(prob) if prob > 0 else -float("inf")

            if log_likelihood > best_likelihood:
                best_likelihood = log_likelihood
                best_config = newconfig

        for h in range(len(haplotypes)):
            haplotypes[h][pos] = best_config[h]

    return haplotypes


def select_clusters(allele_depths, ploidy, max_gap):
    """
    For every position, computes a list of relevant clusters for the threading
    algorithm. Relevant means, that the relative coverage is at least 1/8 of
    what a single haplotype is expected to have for the given ploidy. Apart
    from that, at least <ploidy> and at most <ploidy + 2> many clusters are
    selected to avoid exponential blow-up.
    """
    cov_map = [[] for _ in range(len(allele_depths))]
    for pos in range(len(allele_depths)):
        sorted_cids = sorted(
            ((cid, sum(allele_depths[pos][cid].values())) for cid in allele_depths[pos]),
            key=lambda x: x[1],
            reverse=True,
        )
        total_cov = sum(e[1] for e in sorted_cids)
        cut_off = min(len(sorted_cids), ploidy + 2)
        cov_map[pos].append(sorted_cids[0][0])
        for cid, cov in sorted_cids[1:cut_off]:
            if cov / total_cov < (1.0 / (8.0 * ploidy)) and cov_map[pos]:
                break
            else:
                cov_map[pos].append(cid)

    # re-add clusters missing on at most max_gap intermediate positions
    cut_off = ploidy + 2
    for pos in range(1, len(cov_map) - 1):
        for cid in cov_map[pos - 1]:
            if len(cov_map[pos]) >= cut_off:
                break
            if cid in cov_map[pos]:
                continue
            if any(
                [cid in cov_map[pos + k + 1] for k in range(min(max_gap, len(cov_map) - pos - 1))]
            ):
                cov_map[pos].append(cid)
                allele_depths[pos][cid] = dict()

    for sub in cov_map:
        sub.sort()

    return cov_map


def get_allele_depths(allele_matrix, clustering, ploidy):
    """
    Returns a list, which for every position contains a list (representing the clusters) of dictionaries containing the allele depths.
    Additionally computes a consensus list per position per cluster, such that the first k elements represent the alleles of this
    cluster, if it was selected with multiplicity k.
    ad[pos][c_id][al] = number of reads in cluster c_id having allele al at position pos
    Indices are local, i.e. the i-th entry of ad is the entry for the i-th position that occurs in the allele_matrix.
    """

    num_vars = allele_matrix.getNumPositions()

    # stores allele depth per position and cluster
    ad = [dict() for pos in range(num_vars)]
    cons_lists = [dict() for pos in range(num_vars)]

    # count alleles
    for c_id, cluster in enumerate(clustering):
        for read in cluster:
            for pos, allele in allele_matrix.getRead(read):
                if c_id not in ad[pos]:
                    ad[pos][c_id] = dict()
                if allele not in ad[pos][c_id]:
                    ad[pos][c_id][allele] = 0
                ad[pos][c_id][allele] += 1

    # compute allele lists
    for pos in range(num_vars):
        for c_id in ad[pos]:
            cons_lists[pos][c_id] = []
            cnts = defaultdict(int)
            for i in range(ploidy):
                max_cnt = 0
                max_al = 0
                for al in ad[pos][c_id]:
                    cnt = ad[pos][c_id][al] / (1 + cnts[al])
                    if cnt > max_cnt:
                        max_cnt = cnt
                        max_al = al
                cons_lists[pos][c_id].append(max_al)
                cnts[max_al] += 1

    return ad, cons_lists
