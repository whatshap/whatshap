import logging

from math import log
from collections import defaultdict
from typing import List
from scipy.stats import binom

from whatshap.core import (
    Genotype,
    TriangleSparseMatrix,
)
from whatshap.polyphaseplots import create_histogram
from whatshap.vcf import VariantTable

logger = logging.getLogger(__name__)


def get_variant_scoring(
    variant_table: VariantTable, parent: str, co_parent: str, offspring: List[str], phasing_param
):
    scoring = TriangleSparseMatrix()
    max_dist = 80

    if phasing_param.ploidy != 4:
        return scoring

    ref, alt, alt_count, alt_count_co = classify_variants(variant_table, parent, co_parent)

    num_nodes = 0
    node_to_variant = dict()
    node_positions = []
    simplex_nulliplex_nodes = []
    allowed_pairs = [(1, 0)]
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
            num_nodes += 1

    logger.info("Number of nodes to cluster: %d", num_nodes)
    logger.info("Number of simplex-nulliplex variants: %d", len(simplex_nulliplex_nodes))

    # compute genotype likelihoods for offspring per variant
    scores = []
    off_gl = dict()
    for i, off in enumerate(offspring):
        print("Compute GL for offspring {} out of {}".format(i, len(offspring)))
        off_gl[off], cov = compute_gt_likelihoods(
            variant_table, off, node_positions, ref, alt, 0.01
        )

    for i in range(num_nodes):
        ni = node_to_variant[i]
        # skip if position i is not simplex-nulliplex
        if alt_count[ni] != 1 or alt_count_co[ni] != 0:
            continue

        print(
            "scoring node {}: ref={}, alt={}, count = {} / {}".format(
                i, ref[ni], alt[ni], alt_count[ni], alt_count_co[ni]
            )
        )
        # iterate over next max_dist relevant positions
        for j in range(i + 1, min(i + max_dist + 1, num_nodes)):
            nj = node_to_variant[j]
            if node_to_variant[i] == node_to_variant[j]:
                score = -float("inf")
            else:
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
                    """
                    elif alt_count[nj] == 2 and alt_count_co[nj] == 0:
                        score = score_duplex_nulliplex_tetra(off_gts, off_ad, ref_1=ref[ni], alt_1=alt[ni],  ref_2=ref[nj], alt_2=alt[nj])
                    elif alt_count[nj] == 3 and alt_count_co[nj] == 0:
                        score = score_triplex_nulliplex_tetra(off_gts, off_ad, ref_1=ref[ni], alt_1=alt[ni],  ref_2=ref[nj], alt_2=alt[nj])
                    elif alt_count[nj] == 1 and alt_count_co[nj] == 1:
                        score = score_simplex_simplex_tetra(off_gts, off_ad, ref_1=ref[ni], alt_1=alt[ni],  ref_2=ref[nj], alt_2=alt[nj])
                    """

            scoring.set(i, j, score)
            if score != -float("inf"):
                scores.append(score)

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

    return scoring, node_to_variant


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


def compute_gt_likelihoods(
    variant_table: VariantTable,
    offspring: str,
    positions: List[int],
    ref_alleles: List[int],
    alt_alleles: List[int],
    err: float,
):
    gt_likelihoods = []
    coverages = []
    allele_depths = variant_table.allele_depths_of(offspring)
    for pos in positions:
        gl = dict()
        ref = ref_alleles[pos]
        alt = alt_alleles[pos]
        ref_dp = allele_depths[pos][ref]
        alt_dp = allele_depths[pos][alt]
        coverages.append(ref_dp + alt_dp)
        if ref_dp + alt_dp >= 2:
            gl = defaultdict(float)
            # likelihoods for first variant
            # P(GT=0/0/0/0|allele_depth) = P(allele_depth|GT=0/0/0/0) * P(GT=0/0/0/0) / P(allele_depth)
            # P(GT=0/0/0/1|allele_depth) = P(allele_depth|GT=0/0/0/1) * P(GT=0/0/0/1) / P(allele_depth)
            gl_0 = binom.pmf(alt_dp, ref_dp + alt_dp, err)  # * 0.5
            gl_1 = binom.pmf(alt_dp, ref_dp + alt_dp, 0.75 * err + 0.25 * (1 - err))  # * 0.5
            gl[0] = gl_0  # / (gl_0+gl_1)
            gl[1] = gl_1  # / (gl_0+gl_1)

        gt_likelihoods.append(gl)
    return gt_likelihoods, coverages


def score_simplex_nulliplex_tetra(off_gl, pos1, pos2):
    # initialize with prior probabilities for each phasing
    log_likelihood_cooccur = log(1 / 4)
    log_likelihood_disjoint = log(3 / 4)

    for offspring in off_gl:
        if len(off_gl[offspring][pos1]) > 0 and len(off_gl[offspring][pos2]) > 0:
            likelihood_cooccur = 0.0
            likelihood_disjoint = 0.0

            # P(co-occur|GT) = P(GT|co-occur)*P(co-occur)/P(GT)
            # P(disjoint|GT) = P(GT|disjoint)*P(disjoint)/P(GT)
            for gt, exp in zip(
                [(0, 0), (0, 1), (1, 0), (1, 1)],
                [(1 / 2, 1 / 6), (0, 1 / 3), (0, 1 / 3), (1 / 2, 1 / 6)],
            ):
                gl = off_gl[offspring][pos1][gt[0]] * off_gl[offspring][pos2][gt[1]]
                likelihood_cooccur += gl * exp[0]
                likelihood_disjoint += gl * exp[1]
            log_likelihood_cooccur += log(likelihood_cooccur)
            log_likelihood_disjoint += log(likelihood_disjoint)

    return log_likelihood_cooccur - log_likelihood_disjoint


def score_duplex_nulliplex_tetra(off_gts, off_ad, ref_1=0, alt_1=1, ref_2=0, alt_2=1):
    target_gts = (Genotype([ref_1] * 4), Genotype([ref_2] * 4))
    count_0_0 = count_target_gt(off_gts, target_gts)

    # if both variants share 1-allele: support = 1/6
    pval_e = binom.pmf(count_0_0, len(off_gts), 1 / 6)

    # if both variants don't share 1-allele: support = 0
    pval_d = binom.pmf(count_0_0, len(off_gts), 1 / 24)

    return log(pval_e / pval_d)


def score_triplex_nulliplex_tetra(off_gts, off_ad, ref_1=0, alt_1=1, ref_2=0, alt_2=1):
    target_gts = (Genotype([alt_1] * 4), Genotype([alt_2] * 4))
    count_0_0 = count_target_gt(off_gts, target_gts)

    # if both variants share 1-allele: support = 1/6
    pval_e = binom.pmf(count_0_0, len(off_gts), 1 / 6)

    # if both variants don't share 1-allele: support = 1/2
    pval_d = binom.pmf(count_0_0, len(off_gts), 1 / 2)

    return log(pval_e / pval_d)


def score_simplex_simplex_tetra(off_gts, off_ad, ref_1=0, alt_1=1, ref_2=0, alt_2=1):
    target_gts = (Genotype([ref_1] * 4), Genotype([ref_2] * 4))
    count_0_0 = count_target_gt(off_gts, target_gts)

    # if both variants share 1-allele: support = 1/4
    pval_e = binom.pmf(count_0_0, len(off_gts), 1 / 4)

    # if both variants don't share 1-allele: support = 1/12
    pval_d = binom.pmf(count_0_0, len(off_gts), 1 / 12)

    return log(pval_e / pval_d)


def count_target_gt(off_gts, target_gts):
    count = 0
    for gt1, gt2 in off_gts:
        if gt1 == target_gts[0] and gt2 == target_gts[1]:
            count += 1
    return count
