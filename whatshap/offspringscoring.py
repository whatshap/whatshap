import logging

from math import log
from collections import defaultdict
from typing import List, Iterable, Tuple
from scipy.stats import binom
from scipy.special import binom as binom_coeff

from whatshap.core import TriangleSparseMatrix, ProgenyGenotypeLikelihoods
from whatshap.vcf import VariantTable
from whatshap.variantselection import VariantInfo

logger = logging.getLogger(__name__)


class CachedBinomialCalculator:
    def __init__(self, ploidy, error_rate):
        self.binom_cache = [dict() for _ in range(ploidy + 1)]
        self.ploidy = ploidy
        self.error_rate = error_rate
        self.prob = [
            (1 - g / ploidy) * error_rate + (g / ploidy) * (1 - error_rate)
            for g in range(ploidy + 1)
        ]

    def get_binom_pmf(self, n, k, g):
        if g < 0 or g > self.ploidy or not isinstance(g, int):
            raise ValueError("Invalid genotype alt-count ({}).".format(g))
        if (n, k) not in self.binom_cache[g]:
            self.binom_cache[g][(n, k)] = binom.pmf(k, n, self.prob[g])
        return self.binom_cache[g][(n, k)]


def hyp(N, M, n, k):
    return binom_coeff(M, k) * binom_coeff(N - M, n - k) / binom_coeff(N, n)


def correct_variant_types(
    variant_table: VariantTable,
    progeny_table: VariantTable,
    offspring: List[str],
    varinfo: VariantInfo,
    phasing_param,
):
    # compute unbiased progeny genotype likelihoods
    off_gl = get_offspring_gl(variant_table, progeny_table, offspring, varinfo, phasing_param)
    correction = dict()

    # precompute dosage distributions for every variant type
    k = phasing_param.ploidy
    dosages = []
    for i in range(k + 1):
        for j in range(i + 1):
            d = [
                sum([hyp(k, i, k // 2, l) * hyp(k, j, k // 2, m - l) for l in range(m + 1)])
                for m in range(k + 1)
            ]
            dosages.append(((i, j), tuple(d)))

    var_id = -1
    # compute best fitting variant type based on progeny genotypes
    removing = []
    for node_id in range(off_gl.getNumPositions()):
        if var_id == varinfo.node_to_variant(node_id):
            continue
        var_id = varinfo.node_to_variant(node_id)
        genpos = variant_table.variants[var_id].position
        gt = get_most_likely_variant_type(dosages, genpos, off_gl, node_id)
        alt = varinfo[var_id].alt_count
        co_alt = varinfo[var_id].co_alt_count
        varinfo.correct_type(var_id, gt[0], gt[1])
        if (alt, co_alt) not in correction:
            correction[(alt, co_alt)] = defaultdict(int)
        correction[(alt, co_alt)][gt] += 1
        if not check_variant_compatibility(alt, co_alt, gt[0], gt[1]):
            removing.append(var_id)

    # remove changed variants incompatible to old type (do this after!)
    for var_id in removing:
        varinfo.remove_phasable(var_id)

    # show variant corrections:
    logger.info("   Correcting variant type based on progenies:")
    for old_gt in correction:
        total = sum([correction[old_gt][new_gt] for new_gt in correction[old_gt]])
        if total == 0:
            continue
        logger.info("   {}/{} ({})".format(old_gt[0], old_gt[1], total))
        for new_gt in correction[old_gt]:
            logger.info(
                "      -> {}/{}: {} ({:2.1f}%)".format(
                    new_gt[0],
                    new_gt[1],
                    correction[old_gt][new_gt],
                    100 * correction[old_gt][new_gt] / total,
                )
            )


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
    binom_calc = CachedBinomialCalculator(phasing_param.ploidy, phasing_param.allele_error_rate)
    for i, off in enumerate(offspring):
        gls = compute_gt_likelihoods(
            progeny_table,
            off,
            zip(varinfo.get_node_positions(), progeny_positions),
            varinfo,
            phasing_param.ploidy,
            binom_calc,
            gt_gl_priors,
        )
        for pos, gl in enumerate(gls):
            if gl:
                off_gl.setGlv(pos, i, gl)
    del binom_calc

    return off_gl


def check_variant_compatibility(old_alt, old_co_alt, new_alt, new_co_alt):
    if old_alt == 1 and old_co_alt == 0:
        return (new_alt, new_co_alt) in [(1, 0), (1, 1), (2, 0)]
    elif old_alt == 1 and old_co_alt == 1:
        return (new_alt, new_co_alt) in [(1, 1)]
    elif old_alt == 2 and old_co_alt == 0:
        return (new_alt, new_co_alt) in [(1, 0), (1, 1), (2, 0)]
    return False


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

            scoring.set(i, j, score)

    return scoring


def get_most_likely_variant_type(dosages, genpos, off_gl, pos):

    best_gts = dosages[0][0]
    best_llh = -float("inf")
    for gts, dosage in dosages:
        llh = 1.0
        for i in range(off_gl.getNumSamples()):
            if off_gl.getGl(pos, i, 0) < 0.0:
                continue
            likelihood = 0.0
            for g in range(len(dosage)):
                likelihood += dosage[g] * off_gl.getGl(pos, i, g)
            if likelihood <= 0.0:
                llh -= float("inf")
            else:
                llh += log(likelihood)
        if llh > best_llh:
            best_gts = gts
            best_llh = llh

    return best_gts


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

    prior_dual = [[[0.0] * (ploidy + 1) for _ in range(ploidy + 1)] for _ in range(ploidy + 1)]
    for num_alts_parent in range(0, ploidy + 1):
        for num_alts_coparent in range(0, ploidy + 1):
            for i in range(max_alts + 1):
                for j in range(max_alts + 1):
                    num_alts_offspring = i + j
                    prior_dual[num_alts_parent][num_alts_coparent][num_alts_offspring] += (
                        prior_single[num_alts_parent][i] * prior_single[num_alts_coparent][j]
                    )

    return prior_dual


def compute_gt_likelihoods(
    progeny_table: VariantTable,
    offspring: str,
    position_pairs: Iterable[Tuple[int, int]],
    varinfo: VariantInfo,
    ploidy: int,
    binom_calc: CachedBinomialCalculator,
    gt_priors=None,
):
    gt_likelihoods = []
    allele_depths = progeny_table.allele_depths_of(offspring)

    prev_pos = -1

    for parent_pos, progeny_pos in position_pairs:
        if progeny_pos == prev_pos:
            gt_likelihoods.append(gt_likelihoods[-1])
            continue
        gl = defaultdict(float)
        gl = [0.0 for _ in range(0, ploidy + 1)]
        ref = varinfo[parent_pos].ref
        alt = varinfo[parent_pos].alt
        ref_dp = allele_depths[progeny_pos][ref] if len(allele_depths[progeny_pos]) > ref else 0
        alt_dp = allele_depths[progeny_pos][alt] if len(allele_depths[progeny_pos]) > alt else 0
        num_alts_parent = varinfo[parent_pos].alt_count
        num_alts_coparent = varinfo[parent_pos].co_alt_count
        if ref_dp + alt_dp >= ploidy:
            for i in range(0, ploidy + 1):
                gl[i] = binom_calc.get_binom_pmf(ref_dp + alt_dp, alt_dp, i)
                if gt_priors:
                    gl[i] *= gt_priors[num_alts_parent][num_alts_coparent][i]
            # normalizing likelihoods to sum up to 1 is not necessary, because we compute likelihood
            # ratios later anyways. otherwise it would be done here
            sum_gl = sum(gl)
            for i in range(0, ploidy + 1):
                gl[i] = gl[i] / sum_gl
        else:
            gl = None
        gt_likelihoods.append(gl)
        prev_pos = progeny_pos

    del allele_depths
    return gt_likelihoods
