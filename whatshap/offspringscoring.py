import logging

from collections import defaultdict
from typing import List, Iterable, Tuple
from scipy.stats import binom
from scipy.special import binom as binom_coeff

from whatshap.core import TriangleSparseMatrix, ProgenyGenotypeLikelihoods
from whatshap.vcf import VariantTable, PloidyError

logger = logging.getLogger(__name__)


class CachedBinomialCalculator:

    def __init__(self, ploidy, error_rate):
        self.binom_cache = [dict() for _ in range(ploidy + 1)]
        self.ploidy = ploidy
        self.error_rate = error_rate
        self.prob = [(1 - g / ploidy) * error_rate + (g / ploidy) * (1 - error_rate) for g in range(ploidy + 1)]

    def get_binom_pmf(self, n, k, g):
        if g < 0 or g > self.ploidy or not isinstance(g, int):
            raise ValueError('Invalid genotype alt-count ({}).'.format(g))
        if (n, k) not in self.binom_cache[g]:
            self.binom_cache[g][(n, k)] = binom.pmf(k, n, self.prob[g])
        return self.binom_cache[g][(n, k)]


class VariantInfo:
    __slots__ = ("ref", "alt", "alt_count", "co_alt_count")

    def __init__(self, ref, alt, alt_count, co_alt_count):
        self.ref = ref
        self.alt = alt
        self.alt_count = alt_count
        self.co_alt_count = co_alt_count


def get_phasable_parent_variants(variant_table: VariantTable, parent: str, co_parent: str, phasing_param):
    if phasing_param.ploidy % 2 != 0:
        logger.error("Odd ploidy not supported!")
        raise PloidyError("Odd ploidy not supported!")

    # determine phasable variants
    varinfo = classify_variants(variant_table, parent, co_parent)
    phasable_indices = []
    if phasing_param.simplex:
        allowed_pairs = [(1, 0)]
    else:
        allowed_pairs = [(1, 0), (2, 0), (1, 1)]
    for i in range(len(variant_table)):
        if varinfo[i].alt is None:
            continue
        if (varinfo[i].alt_count, varinfo[i].co_alt_count) not in allowed_pairs:
            continue
        phasable_indices.append(i)

    return varinfo, phasable_indices


def get_offspring_gl(
    variant_table: VariantTable, progeny_table: VariantTable, offspring: List[str], varinfo: List[VariantInfo], phasable_indices: List[int], phasing_param, allele_error_rate=0.06
):

    if phasing_param.ploidy % 2 != 0:
        logger.error("Odd ploidy not supported!")
        raise PloidyError("Odd ploidy not supported!")

    # create map to find genetic positions in progeny table
    genpos_to_progenypos = dict()
    for i in range(len(progeny_table)):
        genpos = progeny_table.variants[i].position
        if genpos:
            genpos_to_progenypos[genpos] = i

    # determine phasable progeny variants
    num_nodes = 0
    node_to_variant = dict()
    type_of_node = []
    node_positions = []
    progeny_positions = []
    simplex_nulliplex_nodes = 0
    for i in phasable_indices:
        genpos = variant_table.variants[i].position
        if genpos not in genpos_to_progenypos:
            continue
        if varinfo[i].alt_count == 1 and varinfo[i].co_alt_count == 0:
            simplex_nulliplex_nodes += 1
        for j in range(varinfo[i].alt_count):
            node_to_variant[num_nodes] = i
            node_positions.append(i)
            progeny_positions.append(genpos_to_progenypos[genpos])
            type_of_node.append((varinfo[i].alt_count, varinfo[i].co_alt_count))
            num_nodes += 1

    logger.info("   Number of nodes to cluster: %d", num_nodes)
    logger.info("   Number of simplex-nulliplex variants: %d", simplex_nulliplex_nodes)

    # compute genotype likelihoods for offspring per variant
    gt_gl_priors = compute_gt_likelihood_priors(phasing_param.ploidy)
    off_gl = ProgenyGenotypeLikelihoods(phasing_param.ploidy, len(offspring), len(node_positions))
    binom_calc = CachedBinomialCalculator(phasing_param.ploidy, allele_error_rate)
    for i, off in enumerate(offspring):
        gls = compute_gt_likelihoods(
            progeny_table,
            off,
            zip(node_positions, progeny_positions),
            varinfo,
            gt_gl_priors,
            phasing_param.ploidy,
            binom_calc,
        )
        for pos, gl in enumerate(gls):
            if gl:
                off_gl.setGlv(pos, i, gl)
    del binom_calc

    return off_gl, node_to_variant, type_of_node


def get_variant_scoring(varinfo, off_gl, node_to_variant, phasing_param):

    num_nodes = len(node_to_variant)
    scoring = TriangleSparseMatrix()
    for i in range(num_nodes):
        ni = node_to_variant[i]

        # iterate over next max_dist relevant positions
        prev_variant = -1
        prev_score = 0
        for j in range(i + 1, min(i + phasing_param.scoring_window + 1, num_nodes)):
            nj = node_to_variant[j]
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
                        score = score_simplex_nulliplex_tetra(off_gl, i, j)
                    elif varinfo[nj].alt_count == 2 and varinfo[nj].co_alt_count == 0:
                        score = score_duplex_nulliplex_tetra(off_gl, i, j)
                    elif varinfo[nj].alt_count == 1 and varinfo[nj].co_alt_count == 1:
                        score = score_simplex_simplex_tetra(off_gl, i, j)
                    prev_score = score
                    prev_variant = nj

            scoring.set(i, j, score)

    return scoring


def classify_variants(variant_table: VariantTable, parent: str, co_parent: str):

    #ref, alt, alt_count, alt_count_co = dict(), dict(), dict(), dict()
    gts1 = variant_table.genotypes_of(parent)
    gts2 = variant_table.genotypes_of(co_parent)
    varinfo = []

    for i in range(len(variant_table)):
        gt1 = gts1[i]
        gt2 = gts2[i]
        gt1v = gt1.as_vector()
        gt2v = gt2.as_vector()

        if gt1.is_none() or gt2.is_none():
            varinfo.append(VariantInfo(None, None, 0, 0))
            continue

        # check homozygosity
        if gt1.is_homozygous():
            varinfo.append(VariantInfo(gt1v[0], None, 0, 0))
            continue

        # count different alleles
        alleles = set()
        for gt in [gt1v, gt2v]:
            for a in gt:
                alleles.add(a)

        alleles = sorted(list(alleles))

        if len(alleles) > 2:
            # genotypes are not bi-allelic
            varinfo.append(VariantInfo(None, None, 0, 0))
            continue

        assert len(alleles) == 2

        # determine majority allele and its multiplicity
        gt1v.sort()
        ref = gt1v[int(len(gt1v) / 2 - 1)]
        alt = gt1v[0] if gt1v[0] != ref else gt1v[-1]
        alt_count = sum([1 if a == alt else 0 for a in gt1v])
        co_alt_count = sum([1 if a == alt else 0 for a in gt2v])
        varinfo.append(VariantInfo(ref, alt, alt_count, co_alt_count))

    return varinfo


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
    position_pairs: Iterable[Tuple[int,int]],
    varinfo: List[VariantInfo],
    gt_priors,
    ploidy: int,
    binom_calc: CachedBinomialCalculator,
):
    gt_likelihoods = []
    allele_depths = progeny_table.allele_depths_of(offspring)

    prev_pos = -1

    for parent_pos, progeny_pos in position_pairs:
        if progeny_pos == prev_pos:
            gt_likelihoods.append(gt_likelihoods[-1])
            continue
        gl = defaultdict(float)
        gl = [0.0 for _ in range(0, ploidy+1)]
        ref = varinfo[parent_pos].ref
        alt = varinfo[parent_pos].alt
        ref_dp = allele_depths[progeny_pos][ref] if len(allele_depths[progeny_pos]) > ref else 0
        alt_dp = allele_depths[progeny_pos][alt] if len(allele_depths[progeny_pos]) > alt else 0
        num_alts_parent = varinfo[parent_pos].alt_count
        num_alts_coparent = varinfo[parent_pos].co_alt_count
        if ref_dp + alt_dp >= ploidy:
            for i in range(0, ploidy + 1):
                gl[i] = binom_calc.get_binom_pmf(ref_dp + alt_dp, alt_dp, i) * gt_priors[num_alts_parent][num_alts_coparent][i]
            # normalizing likelihoods to sum up to 1 is not necessary, because we compute likelihood
            # ratios later anyways. otherwise it would be done here
            #sum_gl = sum(gl)
            #for i in range(0, ploidy + 1):
            #    gl[i] = gl[i] / sum_gl
        else:
            gl = None
        gt_likelihoods.append(gl)
        prev_pos = progeny_pos

    del allele_depths
    return gt_likelihoods


def score_simplex_nulliplex_tetra(off_gl, pos1, pos2):
    return off_gl.getLogLikelihoodDifference(
        pos1, 
        pos2, 
        [(0, 0), (0, 1), (1, 0), (1, 1)], 
        [(1 / 2, 1 / 6), (0, 1 / 3), (0, 1 / 3), (1 / 2, 1 / 6)]
    )


def score_duplex_nulliplex_tetra(off_gl, pos1, pos2):
    return off_gl.getLogLikelihoodDifference(
        pos1, 
        pos2, 
        [(0, 0), (0, 1), (1, 0), (1, 1), (0, 2), (1, 2)], 
        [(1 / 6, 0), (1 / 3, 1 / 3), (0, 1 / 6), (1 / 3, 1 / 3), (0, 1 / 6), (1 / 6, 0)]
    )


def score_simplex_simplex_tetra(off_gl, pos1, pos2):
    return off_gl.getLogLikelihoodDifference(
        pos1,
        pos2,
        [(0, 0), (0, 1), (1, 0), (1, 1), (0, 2), (1, 2)],
        [(1 / 4, 1 / 12), (1 / 4, 1 / 4), (0, 1 / 6), (1 / 4, 1 / 4), (0, 1 / 6), (1 / 4, 1 / 12)],
    )
