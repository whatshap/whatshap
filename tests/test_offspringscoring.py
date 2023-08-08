from scipy.special import binom as binom_coeff

from whatshap.cli.polyphasegenetic import PolyphaseGeneticParameter
from whatshap.polyphase.variantselection import compute_phasable_variants
from whatshap.polyphase.offspringscoring import (
    compute_gt_likelihood_priors,
    compute_gt_likelihoods,
    correct_variant_types,
)
from whatshap.vcf import VcfReader


def old_likelihood_prior_function(ploidy):
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


def test_gt_likelihood_priors():
    for k in range(2, 11):
        priors = compute_gt_likelihood_priors(k)
        priors_old = old_likelihood_prior_function(k)
        for i in range(k + 1):
            for j in range(k + 1):
                for l in range(k + 1):
                    assert abs(priors[i][j][l] - priors_old[i][j][l]) < 0.00000000000001


def test_correct_variant_types():
    table = list(
        VcfReader(
            "tests/data/polyphasegenetic.test.parents.vcf",
            only_snvs=False,
            genotype_likelihoods=False,
            ploidy=4,
            mav=True,
        )
    )[0]
    ptable = list(
        VcfReader(
            "tests/data/polyphasegenetic.test.progeny.vcf.gz",
            only_snvs=False,
            genotype_likelihoods=False,
            ploidy=4,
            mav=True,
            allele_depth=True,
        )
    )[0]

    param = PolyphaseGeneticParameter(4, 20, 0.06, 0, 0, True, True, False, "")

    vi = compute_phasable_variants(table, "Parent_A", "Parent_B", param)
    p1 = set(vi.get_phasable())
    correct_variant_types(table, ptable, ptable.samples, vi, param)
    p2 = vi.get_phasable()
    new_np = [x for x in p1 if x not in p2]
    true_new_np = [18, 21, 30, 35, 37, 51, 69, 71, 98, 107, 110]
    true_new_np += [111, 112, 113, 114, 115, 126, 127, 128]
    assert new_np == true_new_np

    ptable_positions = [v.position for v in ptable.variants]
    for pos in new_np:
        g0 = vi[pos].alt_count
        g1 = vi[pos].co_alt_count
        assert (g0, g1) != (1, 0) or table.variants[pos].position not in ptable_positions


def test_compute_gt_likelihoods():
    table = list(
        VcfReader(
            "tests/data/polyphasegenetic.test.parents.vcf",
            only_snvs=False,
            genotype_likelihoods=False,
            ploidy=4,
            mav=True,
        )
    )[0]
    ptable = list(
        VcfReader(
            "tests/data/polyphasegenetic.test.progeny.vcf.gz",
            only_snvs=False,
            genotype_likelihoods=False,
            ploidy=4,
            mav=True,
            allele_depth=True,
        )
    )[0]

    param = PolyphaseGeneticParameter(4, 20, 0.06, 0, 0, True, True, False, "")
    vi = compute_phasable_variants(table, "Parent_A", "Parent_B", param)
    priors = compute_gt_likelihood_priors(param.ploidy)

    genpos_to_progenypos = dict()
    for i in range(len(ptable)):
        genpos = ptable.variants[i].position
        if genpos:
            genpos_to_progenypos[genpos] = i

    progeny_positions = []
    for i, p in enumerate(vi.get_phasable()):
        genpos = table.variants[p].position
        if genpos not in genpos_to_progenypos:
            vi.remove_phasable(p)

    for p in vi.get_phasable():
        genpos = table.variants[p].position
        alt = vi[p].alt_count
        for j in range(alt):
            progeny_positions.append(genpos_to_progenypos[genpos])

    gls = compute_gt_likelihoods(
        ptable,
        ptable.samples[0],
        zip(vi.get_node_positions(), progeny_positions),
        vi,
        param,
        priors,
    )

    assert gls[0][1] == max(gls[0])
    assert gls[1][1] == max(gls[1])
    assert gls[2][0] == max(gls[2])
    assert gls[3][1] == max(gls[3])
    assert gls[4][1] == max(gls[4])
    assert gls[5][1] == max(gls[5])
    assert gls[6][1] == max(gls[6])
    assert gls[7][1] == max(gls[7])
    assert gls[8][0] == max(gls[8])
    assert gls[15][1] == max(gls[15])
    assert gls[16] is None
    assert gls[17] is None
    assert gls[18][0] == max(gls[18])
    assert gls[-2] == gls[-1]
