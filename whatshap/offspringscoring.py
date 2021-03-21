import logging

from math import log
from collections import defaultdict
from typing import List
from scipy.stats import binom
from scipy.special import binom as binom_coeff

from whatshap.core import TriangleSparseMatrix
from whatshap.polyphaseplots import create_histogram
from whatshap.vcf import VariantTable

logger = logging.getLogger(__name__)


def get_variant_scoring(
    variant_table: VariantTable, parent: str, co_parent: str, offspring: List[str], phasing_param
):
    scoring = TriangleSparseMatrix()
    max_dist = phasing_param.scoring_window

    if phasing_param.ploidy % 2 != 0:
        logger.error("Odd ploidy not supported!")
        return scoring

    ref, alt, alt_count, alt_count_co = classify_variants(variant_table, parent, co_parent)

    num_nodes = 0
    node_to_variant = dict()
    type_of_node = []
    node_positions = []
    simplex_nulliplex_nodes = []
    if phasing_param.simplex:
        allowed_pairs = [(1, 0)]
    else:
        allowed_pairs = [(1, 0), (2, 0), (1, 1)]
    for i in range(len(variant_table)):
        if alt[i] is None:
            continue
        if alt_count[i] == 1 and alt_count_co[i] == 0:
            simplex_nulliplex_nodes.append(num_nodes)
        if (alt_count[i], alt_count_co[i]) not in allowed_pairs:
            continue
        for j in range(alt_count[i]):
            node_to_variant[num_nodes] = i
            node_positions.append(i)
            type_of_node.append((alt_count[i], alt_count_co[i]))
            num_nodes += 1

    logger.info("Number of nodes to cluster: %d", num_nodes)
    logger.info("Number of simplex-nulliplex variants: %d", len(simplex_nulliplex_nodes))

    # compute genotype likelihoods for offspring per variant
    gt_gl_priors = compute_gt_likelihood_priors(phasing_param.ploidy)

    scores = []
    covs = []
    off_gl = dict()
    for i, off in enumerate(offspring):
        print("Compute GL for offspring {} out of {}".format(i, len(offspring)))
        off_gl[off], coverages = compute_gt_likelihoods(
            variant_table,
            off,
            node_positions,
            ref,
            alt,
            alt_count,
            alt_count_co,
            gt_gl_priors,
            phasing_param.ploidy,
            0.06,
        )
        covs.extend(coverages)

    for i in range(num_nodes):
        ni = node_to_variant[i]

        if i % 50 == 0:
            print("scored {}/{} nodes".format(i, num_nodes))
        '''
        print(
            "scoring node {}/{}: ref={}, alt={}, count = {} / {}".format(
                i, num_nodes, ref[ni], alt[ni], alt_count[ni], alt_count_co[ni]
            )
        )
        '''
        # iterate over next max_dist relevant positions
        for j in range(i + 1, min(i + max_dist + 1, num_nodes)):
            nj = node_to_variant[j]
            if ni == nj:
                score = -float("inf")
            else:
                # skip if pair (i, j) if i is not simplex-nulliplex
                if alt_count[ni] != 1 or alt_count_co[ni] != 0:
                    continue

                off_gts = []
                off_ad = []
                for off in offspring:
                    gt = variant_table.genotypes_of(off)
                    ad = variant_table.allele_depths_of(off)
                    if gt[ni].is_none() or gt[nj].is_none():
                        continue
                    else:
                        off_gts.append((gt[ni], gt[nj]))
                        off_ad.append((ad[ni], ad[nj]))

                score = 0.0

                if len(off_gts) > 0 or len(off_ad) > 0:
                    if alt_count[nj] == 1 and alt_count_co[nj] == 0:
                        score = score_simplex_nulliplex_tetra(off_gl, i, j)
                    elif alt_count[nj] == 2 and alt_count_co[nj] == 0:
                        score = score_duplex_nulliplex_tetra(off_gl, i, j)
                    elif alt_count[nj] == 1 and alt_count_co[nj] == 1:
                        score = score_simplex_simplex_tetra(off_gl, i, j)

            scoring.set(i, j, score)
            if score != -float("inf"):
                scores.append(score)

    """
    create_histogram(
        "scores{}.pdf".format(len(variant_table)),
        scores,
        [],
        100,
        [min(scores), max(scores)],
        "scores",
        "Scores among Simplex-Nulliplex variant pairs",
        name1="score",
        name2="n/a",
    )

    create_histogram(
        "coverages{}.pdf".format(len(variant_table)),
        covs,
        [],
        100,
        [0, 100],
        "coverages",
        "Offspring coverages among Simplex-Nulliplex variant pairs",
        name1="coverage",
        name2="n/a",
    )
    """

    return scoring, node_positions, type_of_node


def classify_variants(variant_table: VariantTable, parent: str, co_parent: str):

    ref, alt, alt_count, alt_count_co = dict(), dict(), dict(), dict()
    gts1 = variant_table.genotypes_of(parent)
    gts2 = variant_table.genotypes_of(co_parent)

    for i in range(len(variant_table)):
        gt1 = gts1[i]
        gt2 = gts2[i]
        gt1v = gt1.as_vector()
        gt2v = gt2.as_vector()

        if gt1.is_none() or gt2.is_none():
            ref[i] = None
            alt[i] = None
            alt_count[i] = 0
            alt_count_co[i] = 0
            continue

        # check homozygosity
        if gt1.is_homozygous():
            ref[i] = gt1v[0]
            alt[i] = None
            alt_count[i] = 0
            alt_count_co[i] = 0
            continue

        # count different alleles
        alleles = set()
        for gt in [gt1v, gt2v]:
            for a in gt:
                alleles.add(a)

        alleles = sorted(list(alleles))

        if len(alleles) > 2:
            # genotypes are not bi-allelic
            ref[i] = None
            alt[i] = None
            alt_count[i] = 0
            alt_count_co[i] = 0
            continue

        assert len(alleles) == 2

        # determine majority allele and its multiplicity
        gt1v.sort()
        ref[i] = gt1v[int(len(gt1v) / 2 - 1)]
        alt[i] = gt1v[0] if gt1v[0] != ref[i] else gt1v[-1]
        alt_count[i] = sum([1 if a == alt[i] else 0 for a in gt1v])
        alt_count_co[i] = sum([1 if a == alt[i] else 0 for a in gt2v])

    return ref, alt, alt_count, alt_count_co


def compute_gt_likelihood_priors(ploidy):
    # auxiliary table for prior probailities
    max_alts = ploidy // 2  # max. alleles inherited from one parent
    prior_single = [[0.0] * (max_alts + 1) for _ in range(ploidy + 1)]
    for num_alts in range(0, ploidy + 1):
        for num_drawn_alts in range(0, max_alts + 1):
            if ploidy - num_alts >= max_alts - num_drawn_alts and num_alts >= num_drawn_alts:
                prior_single[num_alts][num_drawn_alts] = (
                    binom_coeff(ploidy - num_alts, max_alts - num_drawn_alts)
                    * binom_coeff(num_alts, num_drawn_alts)
                    / binom_coeff(ploidy, max_alts)
                )
            '''
            print(
                "prior[{}][{}] = {}".format(
                    num_alts, num_drawn_alts, prior_single[num_alts][num_drawn_alts]
                )
            )
            '''

    prior_dual = [[[0.0] * (ploidy + 1) for _ in range(ploidy + 1)] for _ in range(ploidy + 1)]
    for num_alts_parent in range(0, ploidy + 1):
        for num_alts_coparent in range(0, ploidy + 1):
            for i in range(max_alts + 1):
                for j in range(max_alts + 1):
                    num_alts_offspring = i + j
                    prior_dual[num_alts_parent][num_alts_coparent][num_alts_offspring] += (
                        prior_single[num_alts_parent][i] * prior_single[num_alts_coparent][j]
                    )
            '''
            for num_alts_offspring in range(0, ploidy + 1):
                print(
                    "prior[{},{}][{}] = {}".format(
                        num_alts_parent,
                        num_alts_coparent,
                        num_alts_offspring,
                        prior_dual[num_alts_parent][num_alts_coparent][num_alts_offspring],
                    )
                )
            '''
    return prior_dual


def compute_gt_likelihoods(
    variant_table: VariantTable,
    offspring: str,
    positions: List[int],
    ref_alleles: List[int],
    alt_alleles: List[int],
    alt_counts: List[int],
    alt_co_counts: List[int],
    gt_priors,
    ploidy: int,
    err: float,
):
    gt_likelihoods = []
    coverages = []
    allele_depths = variant_table.allele_depths_of(offspring)
    for pos in positions:
        gl = defaultdict(float)
        ref = ref_alleles[pos]
        alt = alt_alleles[pos]
        ref_dp = allele_depths[pos][ref]
        alt_dp = allele_depths[pos][alt]
        num_alts_parent = alt_counts[pos]
        num_alts_coparent = alt_co_counts[pos]
        coverages.append(ref_dp + alt_dp)
        if ref_dp + alt_dp >= ploidy:
            gl = defaultdict(float)
            for i in range(0, ploidy + 1):
                gl[i] = (
                    binom.pmf(
                        alt_dp, ref_dp + alt_dp, (1 - i / ploidy) * err + (i / ploidy) * (1 - err)
                    )
                    * gt_priors[num_alts_parent][num_alts_coparent][i]
                )

        gt_likelihoods.append(gl)
    return gt_likelihoods, coverages


def score_simplex_nulliplex_tetra(off_gl, pos1, pos2):
    return score_variant_pair_tetra(
        off_gl,
        pos1,
        pos2,
        [(0, 0), (0, 1), (1, 0), (1, 1)],
        [(1 / 2, 1 / 6), (0, 1 / 3), (0, 1 / 3), (1 / 2, 1 / 6)],
    )


def score_duplex_nulliplex_tetra(off_gl, pos1, pos2):
    return score_variant_pair_tetra(
        off_gl,
        pos1,
        pos2,
        [(0, 0), (0, 1), (1, 0), (1, 1), (0, 2), (1, 2)],
        [(1 / 6, 0), (1 / 3, 1 / 3), (0, 1 / 6), (1 / 3, 1 / 3), (0, 1 / 6), (1 / 6, 0)],
    )


def score_simplex_simplex_tetra(off_gl, pos1, pos2):
    return score_variant_pair_tetra(
        off_gl,
        pos1,
        pos2,
        [(0, 0), (0, 1), (1, 0), (1, 1), (0, 2), (1, 2)],
        [(1 / 4, 1 / 12), (1 / 4, 1 / 4), (0, 1 / 6), (1 / 4, 1 / 4), (0, 1 / 6), (1 / 4, 1 / 12)],
    )


def score_variant_pair_tetra(off_gl, pos1, pos2, gt_combinations, gt_frequencies):
    # initialize with prior probabilities for each phasing
    log_likelihood_cooccur = log(1 / 4)
    log_likelihood_disjoint = log(3 / 4)

    for offspring in off_gl:
        if len(off_gl[offspring][pos1]) > 0 and len(off_gl[offspring][pos2]) > 0:
            likelihood_cooccur = 0.0
            likelihood_disjoint = 0.0

            # P(co-occur|GT) = P(GT|co-occur)*P(co-occur)/P(GT)
            # P(disjoint|GT) = P(GT|disjoint)*P(disjoint)/P(GT)
            for gt, exp in zip(gt_combinations, gt_frequencies):
                gl = off_gl[offspring][pos1][gt[0]] * off_gl[offspring][pos2][gt[1]]
                likelihood_cooccur += gl * exp[0]
                likelihood_disjoint += gl * exp[1]
            log_likelihood_cooccur += log(likelihood_cooccur)
            log_likelihood_disjoint += log(likelihood_disjoint)

    return log_likelihood_cooccur - log_likelihood_disjoint
