import logging

from math import log
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
    __slots__ = ("ref", "alt", "alt_count", "co_alt_count", "alt_count_corrected", "co_alt_count_corrected")

    def __init__(self, ref, alt, alt_count, co_alt_count, alt_count_corrected=None, co_alt_count_corrected=None):
        self.ref = ref
        self.alt = alt
        self.alt_count = alt_count
        self.co_alt_count = co_alt_count
        self.alt_count_corrected = alt_count_corrected
        self.co_alt_count_corrected = co_alt_count_corrected


def get_phasable_parent_variants(variant_table: VariantTable, parent: str, co_parent: str, phasing_param):
    if phasing_param.ploidy % 2 != 0:
        logger.error("Odd ploidy not supported!")
        raise PloidyError("Odd ploidy not supported!")

    # determine phasable variants
    varinfo = classify_variants(variant_table, parent, co_parent)
    
    # determine unbiased progeny genotype likelihoods
    
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


def add_corrected_variant_types(variant_table: VariantTable, progeny_table: VariantTable, offspring: List[str], varinfo: List[VariantInfo], phasable_indices: List[int], phasing_param, allele_error_rate=0.06):
    # compute unbiased progeny genotype likelihoods
    off_gl, node_to_variant, type_of_node = get_offspring_gl(variant_table, progeny_table, offspring, varinfo, phasable_indices, phasing_param, allele_error_rate, biased_gl=True)
    correction = dict()
    
    # compute best fitting variant type based on progeny genotypes
    for node_id in range(off_gl.getNumPositions()):
        var_id = node_to_variant[node_id]
        genpos = variant_table.variants[var_id].position
        gt = get_most_likely_variant_type(genpos, varinfo[var_id], off_gl, node_id)
        varinfo[var_id].alt_count_corrected = gt[0]
        varinfo[var_id].co_alt_count_corrected = gt[1]
        if (varinfo[var_id].alt_count, varinfo[var_id].co_alt_count) not in correction:
            correction[(varinfo[var_id].alt_count, varinfo[var_id].co_alt_count)] = defaultdict(int)
        correction[(varinfo[var_id].alt_count, varinfo[var_id].co_alt_count)][gt] += 1
        
    # show variant corrections:
    print("      Variant type corrections:")
    for old_gt in correction:
        total = sum([correction[old_gt][new_gt] for new_gt in correction[old_gt]])
        if total == 0:
            continue
        print("      {}/{} ({})".format(old_gt[0], old_gt[1], total))
        for new_gt in correction[old_gt]:
            print("         -> {}/{}: {} ({:2.1f}%)".format(new_gt[0], new_gt[1], correction[old_gt][new_gt], 100*correction[old_gt][new_gt]/total))
    
    return varinfo

def get_offspring_gl(
    variant_table: VariantTable, progeny_table: VariantTable, offspring: List[str], varinfo: List[VariantInfo], phasable_indices: List[int], phasing_param, allele_error_rate=0.06, biased_gl=True
):

    if phasing_param.ploidy != 4:
        logger.error("Only ploidy 4 is supported!")
        raise PloidyError("Only ploidy 4 is supported!")

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
        if varinfo[i].alt_count_corrected:
            alt = varinfo[i].alt_count_corrected
            co_alt = varinfo[i].co_alt_count_corrected
            if not check_variant_compatibility(varinfo[i].alt_count, varinfo[i].co_alt_count, alt, co_alt):
                continue
        else:
            alt = varinfo[i].alt_count
            co_alt = varinfo[i].co_alt_count
        if alt == 1 and co_alt == 0:
            simplex_nulliplex_nodes += 1
        for j in range(alt):
            node_to_variant[num_nodes] = i
            node_positions.append(i)
            progeny_positions.append(genpos_to_progenypos[genpos])
            type_of_node.append((alt, co_alt))
            num_nodes += 1

    logger.info("   Number of nodes to cluster: %d", num_nodes)
    logger.info("   Number of simplex-nulliplex variants: %d", simplex_nulliplex_nodes)

    # compute genotype likelihoods for offspring per variant
    if biased_gl:
        gt_gl_priors = compute_gt_likelihood_priors(phasing_param.ploidy)
    else:
        gt_gl_priors = None
    off_gl = ProgenyGenotypeLikelihoods(phasing_param.ploidy, len(offspring), len(node_positions))
    binom_calc = CachedBinomialCalculator(phasing_param.ploidy, allele_error_rate)
    for i, off in enumerate(offspring):
        gls = compute_gt_likelihoods(
            progeny_table,
            off,
            zip(node_positions, progeny_positions),
            varinfo,
            phasing_param.ploidy,
            binom_calc,
            gt_gl_priors,
            use_corrected_counts = not (varinfo[0].alt_count_corrected is None)
        )
        for pos, gl in enumerate(gls):
            if gl:
                off_gl.setGlv(pos, i, gl)
    del binom_calc

    return off_gl, node_to_variant, type_of_node


def check_variant_compatibility(old_alt, old_co_alt, new_alt, new_co_alt):
    if old_alt == 1 and old_co_alt == 0:
        return (new_alt, new_co_alt) in [(1,0), (1,1), (2,0)]
    elif old_alt == 1 and old_co_alt == 1:
        return (new_alt, new_co_alt) in [(1,1)]
    elif old_alt == 2 and old_co_alt == 0:
        return (new_alt, new_co_alt) in [(1,0), (1,1), (2,0)]
    return false

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


def get_most_likely_variant_type(genpos, varinfo, off_gl, pos):
    '''
    if (genpos+1) not in [56711005, 56711020, 56711024, 56711025, 56711037, 56711038, 56711071, 
                          56711072, 56711076, 56711078, 56711083, 56711084, 56711097, 56711110,
                          56711132, 56711161, 56711163, 56711164, 56711171, 56711173, 56711175,
                          56711176, 56711177, 56711179, 56711182, 56711183, 56711185, 56711186, 
                          56711188, 56711189, 56711190, 56711191, 56711192, 56711195, 56711199, 
                          56711211, 56711213, 56711218, 56711230, 56711232, 56711245, 56711274,
                          56711279, 56711290, 56711304, 56711306, 56711322, 56711326, 56711332, 
                          56711333, 56711336, 56711339, 56711340, 56711342, 56711346, 56711347,
                          56711348, 56711349, 56711350, 56711352, 56711353, 56711356, 56711358,
                          56711359, 56711361, 56711367, 56711368, 56711369, 56711371, 56711373,
                          56711374, 56711375, 56711376, 56711378, 56711379, 56711381, 56711383,
                          56711385, 56711386, 56711387, 56711391, 56711392, 56711393, 56711396,
                          56711397, 56711398, 56711399, 56711400, 56711402, 56711403, 56711407,
                          56713030, 56713053, 56713072, 56713079, 56713126, 56713149, 56713150,
                          56713155, 56713196, 56713204, 56713207, 56713245, 56713263, 56713265,
                          56713273, 56713284, 56713303, 56713318, 56713349, 56713351, 56713418,
                          56713428, 56713460, 56713461, 56713467, 56713468, 56713487, 56713494,
                          56713498, 56713505, 56713518, 56713530, 56713538, 56713543, 56713545
                         ]:
        return None
    '''
    
    dosages = [((0,0), (1.0, 0, 0, 0, 0)),
               ((1,0), (0.5, 0.5, 0, 0, 0)),
               ((1,1), (0.25, 0.5, 0.25, 0, 0)),
               ((2,0), (1/6, 4/6, 1/6, 0, 0)),
               ((2,1), (1/12, 5/12, 5/12, 1/12, 0)),
               ((2,2), (1/36, 2/9, 0.5, 2/9, 1/36)),
               ((3,0), (0, 0.5, 0.5, 0, 0)),
               ((3,1), (0, 0.25, 0.5, 0.25, 0)),
               ((3,2), (0, 1/12, 5/12, 5/12, 1/12)),
               ((3,3), (0, 0, 0.25, 0.5, 0.25)),
               ((4,0), (0, 0, 1.0, 0, 0)),
               ((4,1), (0, 0, 0.5, 0.5, 0)),
               ((4,2), (0, 0, 1/6, 4/6, 1/6)),
               ((4,3), (0, 0, 0, 0.5, 0.5)),
               ((4,4), (0, 0, 0, 0, 1.0))]
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
    position_pairs: Iterable[Tuple[int,int]],
    varinfo: List[VariantInfo],
    ploidy: int,
    binom_calc: CachedBinomialCalculator,
    gt_priors=None,
    use_corrected_counts=False
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
        if varinfo[parent_pos].alt_count_corrected and varinfo[parent_pos].co_alt_count_corrected:
            num_alts_parent = varinfo[parent_pos].alt_count_corrected
            num_alts_coparent = varinfo[parent_pos].co_alt_count_corrected
        else:
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
