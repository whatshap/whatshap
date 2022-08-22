from whatshap.core import Pedigree, PhredGenotypeLikelihoods, NumericSampleIds

from whatshap.pedigree import find_recombination, RecombinationEvent
from whatshap.testhelpers import canonic_index_list_to_biallelic_gt_list


def test_pedigree_no_gls():
    ped = Pedigree(NumericSampleIds())
    genotypes1 = canonic_index_list_to_biallelic_gt_list([0, 1, 0, 2])
    genotypes5 = canonic_index_list_to_biallelic_gt_list([1, 2, 2, 0])
    ped.add_individual("sample1", genotypes1)
    assert len(ped) == 1
    assert ped.variant_count == 4, str(ped.variant_count)
    ped.add_individual("sample5", genotypes5)
    assert len(ped) == 2
    assert ped.variant_count == 4, str(ped.variant_count)
    for i in range(ped.variant_count):
        assert ped.genotype("sample1", i) == genotypes1[i]
        assert ped.genotype_likelihoods("sample1", i) is None
        assert ped.genotype("sample5", i) == genotypes5[i]
        assert ped.genotype_likelihoods("sample5", i) is None


def test_pedigree_with_gls():
    ped = Pedigree(NumericSampleIds())
    genotypes1 = canonic_index_list_to_biallelic_gt_list([0, 1, 0, 2])
    gls1 = [
        PhredGenotypeLikelihoods([0, 1, 2]),
        PhredGenotypeLikelihoods([215, 81, 147]),
        PhredGenotypeLikelihoods([199, 49, 253]),
        PhredGenotypeLikelihoods([167, 200, 163]),
    ]
    genotypes5 = canonic_index_list_to_biallelic_gt_list([1, 2, 2, 0])
    gls5 = [
        PhredGenotypeLikelihoods([184, 71, 233]),
        PhredGenotypeLikelihoods([65, 32, 87]),
        PhredGenotypeLikelihoods([28, 215, 131]),
        PhredGenotypeLikelihoods([98, 250, 137]),
    ]
    ped.add_individual("sample1", genotypes1, gls1)
    assert len(ped) == 1
    assert ped.variant_count == 4, str(ped.variant_count)
    ped.add_individual("sample5", genotypes5, gls5)
    assert len(ped) == 2
    assert ped.variant_count == 4, str(ped.variant_count)
    for i in range(ped.variant_count):
        assert ped.genotype("sample1", i) == genotypes1[i]
        assert list(ped.genotype_likelihoods("sample1", i)) == list(gls1[i])
        assert ped.genotype("sample5", i) == genotypes5[i]
        assert list(ped.genotype_likelihoods("sample5", i)) == list(gls5[i])


def test_find_recombination():
    transmission_vector = [0, 0, 1, 1, 0]
    positions = [5303, 5432, 8307, 9000, 9500]
    recombcost = [0, 3, 3, 1, 1]
    components = {5303: 5303, 5432: 5303, 8307: 5303, 9000: 5303, 9500: 5303}
    events = find_recombination(transmission_vector, components, positions, recombcost)
    assert events == [
        RecombinationEvent(
            position1=5432,
            position2=8307,
            transmitted_hap_father1=0,
            transmitted_hap_father2=1,
            transmitted_hap_mother1=0,
            transmitted_hap_mother2=0,
            recombination_cost=3,
        ),
        RecombinationEvent(
            position1=9000,
            position2=9500,
            transmitted_hap_father1=1,
            transmitted_hap_father2=0,
            transmitted_hap_mother1=0,
            transmitted_hap_mother2=0,
            recombination_cost=1,
        ),
    ]
