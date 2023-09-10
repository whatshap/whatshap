"""
Computes allele co-occurence scores for genetic polyploid phasing.

For a set of phasable variants the genetic polyphase algorithm computes scores to quantify how
likely the marker alleles (= minority alleles) of two variants origin from the same of from
different haplotypes. The computation is done in several steps, as described in the
`polyphasegenetic` module.
"""

import logging

from math import log, isnan
from collections import defaultdict
from typing import List, Iterable, Tuple
from scipy.stats import binom
from scipy.special import binom as binom_coeff
from functools import lru_cache

from whatshap.polyphase.solver import TriangleSparseMatrix, ProgenyGenotypeLikelihoods
from whatshap.polyphase.variantselection import VariantInfo
from whatshap.vcf import VariantTable

logger = logging.getLogger(__name__)


@lru_cache(maxsize=None)
def get_binom_pmf(n, k, g, ploidy, error_rate):
    if g < 0 or g > ploidy or not isinstance(g, int):
        raise ValueError(f"Invalid genotype alt-count ({g}).")
    return binom.pmf(k, n, (1 - g / ploidy) * error_rate + (g / ploidy) * (1 - error_rate))


def hyp(k, N, M, n):
    return binom_coeff(M, k) * binom_coeff(N - M, n - k) / binom_coeff(N, n)


def correct_variant_types(
    variant_table: VariantTable,
    progeny_table: VariantTable,
    offspring: List[str],
    varinfo: VariantInfo,
    phasing_param,
):
    # compute unbiased progeny genotype likelihoods
    priors = compute_gt_likelihood_priors(phasing_param.ploidy)
    off_gl = get_offspring_gl(variant_table, progeny_table, offspring, varinfo, phasing_param)
    correction = dict()

    # compute best fitting variant type based on progeny genotypes
    var_id = -1
    correcting = []
    for node_id in range(off_gl.getNumPositions()):
        if var_id == varinfo.node_to_variant(node_id):
            continue

        # determine best fit
        var_id = varinfo.node_to_variant(node_id)
        genpos = variant_table.variants[var_id].position
        gt = get_most_likely_variant_type(priors, genpos, off_gl, node_id)
        correcting.append((var_id, gt))

        # count for statistics
        alt = varinfo[var_id].alt_count
        co_alt = varinfo[var_id].co_alt_count
        if (alt, co_alt) not in correction:
            correction[(alt, co_alt)] = defaultdict(int)
        correction[(alt, co_alt)][gt] += 1

    # apply changes afterwards (do not inside the loop above!)
    for var_id, gt in correcting:
        varinfo.correct_type(var_id, gt[0], gt[1])

    # show variant corrections:
    logger.info("   Correcting variant type based on progenies:")
    for old_gt in correction:
        total = sum([correction[old_gt][new_gt] for new_gt in correction[old_gt]])
        if total == 0:
            continue
        logger.info(f"   {old_gt[0]}/{old_gt[1]} ({total})")
        for new_gt in correction[old_gt]:
            num = correction[old_gt][new_gt]
            perc = 100 * correction[old_gt][new_gt] / total
            logger.info("%s", f"      -> {new_gt[0]}/{new_gt[1]}: {num} ({perc:2.1f}%)")


def get_offspring_gl(
    variant_table: VariantTable,
    progeny_table: VariantTable,
    offspring: List[str],
    varinfo: VariantInfo,
    phasing_param,
):
    # create map to find genetic positions in progeny table
    genpos_to_progenypos = dict()
    for i in range(len(progeny_table)):
        genpos = progeny_table.variants[i].position
        if genpos:
            genpos_to_progenypos[genpos] = i

    # determine phasable progeny variants
    num_nodes = 0
    progeny_positions = []
    simplex_nulliplex_nodes = 0
    for i, p in enumerate(varinfo.get_phasable()):
        genpos = variant_table.variants[p].position
        if genpos not in genpos_to_progenypos:
            varinfo.remove_phasable(p)

    for p in varinfo.get_phasable():
        genpos = variant_table.variants[p].position
        alt = varinfo[p].alt_count
        co_alt = varinfo[p].co_alt_count
        if alt == 1 and co_alt == 0:
            simplex_nulliplex_nodes += 1
        for j in range(alt):
            progeny_positions.append(genpos_to_progenypos[genpos])
            num_nodes += 1

    logger.info("   Number of nodes to cluster: %d", num_nodes)
    logger.info("   Number of simplex-nulliplex variants: %d", simplex_nulliplex_nodes)

    # compute genotype likelihoods for offspring per variant
    gt_gl_priors = compute_gt_likelihood_priors(phasing_param.ploidy)
    off_gl = ProgenyGenotypeLikelihoods(
        phasing_param.ploidy, len(offspring), len(varinfo.get_node_positions())
    )
    for i, off in enumerate(offspring):
        gls = compute_gt_likelihoods(
            progeny_table,
            off,
            zip(varinfo.get_node_positions(), progeny_positions),
            varinfo,
            phasing_param,
            gt_gl_priors,
        )
        for pos, gl in enumerate(gls):
            if gl:
                off_gl.setGlv(pos, i, gl)

    return off_gl


def get_variant_scoring(varinfo, off_gl, phasing_param):
    num_nodes = len(varinfo.get_node_positions())
    scoring = TriangleSparseMatrix()

    # create a stride pattern for the scoring:
    # 25% of samples are direct neighbours, 25% neighbours with stride 3, 7 and 13 each
    w = phasing_param.scoring_window
    w3, w7, w13 = w // 4, w // 2, 3 * w // 4
    strides = [i for i in range(1, w3 + 1)]
    strides += [strides[-1] + 3 * i for i in range(1, w7 - w3 + 1)]
    strides += [strides[-1] + 7 * i for i in range(1, w13 - w7 + 1)]
    strides += [strides[-1] + 13 * i for i in range(1, w - w13 + 1)]

    for i in range(num_nodes):
        ni = varinfo.node_to_variant(i)

        # iterate over next max_dist relevant positions
        prev_variant = -1
        prev_score = 0
        for j in [i + s for s in strides if i + s < num_nodes]:
            nj = varinfo.node_to_variant(j)
            if ni == nj:
                score = -float("inf")
            else:
                # skip if pair (i, j) if i is not simplex-nulliplex
                if varinfo[ni].alt_count != 1 or varinfo[ni].co_alt_count != 0:
                    continue
                if nj == prev_variant:
                    # if j has same position as previous node, if must be multiplex-nulliplex
                    # then the score is identical to the previous node
                    score = prev_score
                else:
                    if varinfo[nj].alt_count == 1 and varinfo[nj].co_alt_count == 0:
                        score = off_gl.getSimplexNulliplexScore(i, j)
                    elif varinfo[nj].alt_count == 2 and varinfo[nj].co_alt_count == 0:
                        score = off_gl.getDuplexNulliplexScore(i, j)
                    elif varinfo[nj].alt_count == 1 and varinfo[nj].co_alt_count == 1:
                        score = off_gl.getSimplexSimplexScore(i, j)
                    prev_score = score
                    prev_variant = nj

            assert score != float("inf")
            assert not isnan(score)
            scoring.set(i, j, score)

    return scoring


def get_most_likely_variant_type(priors, genpos, off_gl, pos):
    best_gts = (0, 0)
    best_llh = -float("inf")
    k = len(priors)
    for g0 in range(k):
        for g1 in range(g0 + 1):
            llh = 1.0
            for i in range(off_gl.getNumSamples()):
                if off_gl.getGl(pos, i, 0) < 0.0:
                    continue
                likelihood = 0.0
                for g in range(k):
                    likelihood += priors[g0][g1][g] * off_gl.getGl(pos, i, g)
                if likelihood <= 0.0:
                    llh -= float("inf")
                else:
                    llh += log(likelihood)
            if llh > best_llh:
                best_gts = (g0, g1)
                best_llh = llh
    return best_gts


def compute_gt_likelihood_priors(ploidy):
    # Semantic: priors[i][j][l] = Prior probability that a progeny inherits l alternative alleles
    # from the parents under the assumption that first parent has i out of k alt. alleles and
    # the second parent has j out of k.
    k = ploidy
    priors = [[[] for _ in range(k + 1)] for _ in range(k + 1)]
    for i in range(k + 1):
        for j in range(i + 1):
            d = [
                sum([hyp(l, k, i, k // 2) * hyp(m - l, k, j, k // 2) for l in range(m + 1)])
                for m in range(k + 1)
            ]
            priors[i][j] = d
            priors[j][i] = d

    return priors


def compute_gt_likelihoods(
    progeny_table: VariantTable,
    offspring: str,
    position_pairs: Iterable[Tuple[int, int]],
    varinfo: VariantInfo,
    param,
    gt_priors=None,
):
    gt_likelihoods = []
    allele_depths = progeny_table.allele_depths_of(offspring)

    prev_pos = -1

    for parent_pos, progeny_pos in position_pairs:
        if progeny_pos == prev_pos:
            gt_likelihoods.append(gt_likelihoods[-1])
            continue
        gl = [0.0 for _ in range(0, param.ploidy + 1)]
        ref = varinfo[parent_pos].ref
        alt = varinfo[parent_pos].alt
        ref_dp = allele_depths[progeny_pos][ref] if len(allele_depths[progeny_pos]) > ref else 0
        alt_dp = allele_depths[progeny_pos][alt] if len(allele_depths[progeny_pos]) > alt else 0
        num_alts_parent = varinfo[parent_pos].alt_count
        num_alts_coparent = varinfo[parent_pos].co_alt_count
        if ref_dp + alt_dp >= param.ploidy:
            for i in range(0, param.ploidy + 1):
                gl[i] = get_binom_pmf(
                    ref_dp + alt_dp, alt_dp, i, param.ploidy, param.allele_error_rate
                )
                if gt_priors:
                    gl[i] *= gt_priors[num_alts_parent][num_alts_coparent][i]
            # normalizing likelihoods to sum up to 1 is not necessary, because we compute likelihood
            # ratios later anyways. otherwise it would be done here
            sum_gl = sum(gl)
            for i in range(0, param.ploidy + 1):
                gl[i] = gl[i] / sum_gl
        else:
            gl = None
        gt_likelihoods.append(gl)
        prev_pos = progeny_pos

    del allele_depths
    return gt_likelihoods
