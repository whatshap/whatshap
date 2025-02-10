import math
from whatshap.core import (
    ReadSet,
    Pedigree,
    NumericSampleIds,
    PhredGenotypeLikelihoods,
    GenotypeDPTable,
)
from whatshap.testhelpers import (
    string_to_readset,
    canonic_index_to_biallelic_gt,
    canonic_index_list_to_biallelic_gt_list,
    likelihoods_equal,
)


def compare_to_expected(dp_forward_backward, positions, expected=None, genotypes=None):
    # check if computed likelihoods are equal to expected ones (if given)
    if expected is not None:
        for i in range(len(positions)):
            likelihoods = dp_forward_backward.get_genotype_likelihoods("individual0", i)
            print(f"Position {i}:")
            print(f"Computed likelihoods: {likelihoods}")
            print(f"Expected likelihoods: {expected[i]}")
            assert likelihoods_equal(
                likelihoods, expected[i]
            ), f"Likelihood mismatch at position {i}: Expected {expected[i]} but got {likelihoods}"

    # check if likeliest genotype is equal to expected genotype
    for i in range(len(positions)):
        likelihoods = dp_forward_backward.get_genotype_likelihoods("individual0", i)
        # find likeliest genotype
        max_val = -1
        max_geno = -1
        for genotype in likelihoods.genotypes():
            # likelihood should not be nan
            assert not math.isnan(likelihoods[genotype])
            if likelihoods[genotype] > max_val:
                max_val = likelihoods[genotype]
                max_geno = genotype

        print(
            "genotype likelihoods for position", i, likelihoods, " likeliest genotype: ", max_geno
        )

        if genotypes is not None:
            assert (
                max_geno == genotypes[i]
            ), f"Mismatch at position {i}: {max_geno} != {genotypes[i]}"


def test_genotyping_empty_readset():
    rs = ReadSet()
    genotypes = canonic_index_list_to_biallelic_gt_list([1, 1])
    recombcost = [1, 1]
    numeric_sample_ids = NumericSampleIds()
    pedigree = Pedigree(numeric_sample_ids)
    genotype_likelihoods = [None, None]
    pedigree.add_individual("individual0", genotypes, genotype_likelihoods)
    _ = GenotypeDPTable(numeric_sample_ids, rs, recombcost, pedigree)


def check_genotyping_single_individual(
    reads, weights=None, expected=None, genotypes=None, scaling=None, genotype_priors=None
):
    # 0) set up read set
    readset = string_to_readset(s=reads, w=weights, scale_quality=scaling)
    positions = readset.get_positions()

    # 1) Genotype using forward backward algorithm
    recombcost = [1] * len(positions)
    numeric_sample_ids = NumericSampleIds()
    pedigree = Pedigree(numeric_sample_ids)
    genotype_likelihoods = [PhredGenotypeLikelihoods([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0])] * len(
        positions
    )

    if genotype_priors is not None:
        genotype_likelihoods = genotype_priors

    pedigree.add_individual(
        "individual0",
        [canonic_index_to_biallelic_gt(1) for i in range(len(positions))],
        genotype_likelihoods,
    )
    dp_forward_backward = GenotypeDPTable(numeric_sample_ids, readset, recombcost, pedigree)

    # check the results
    compare_to_expected(dp_forward_backward, positions, expected, genotypes)


# TODO: when using non-uniform transitions, the pedigree results are not equal (?)
# 2) Phase using PedMEC code for trios with two "empty" individuals (i.e. having no reads)
#    recombcost = [1] * len(positions) # recombination costs 1, should not occur
#    numeric_sample_ids = NumericSampleIds()
#    pedigree = Pedigree(numeric_sample_ids)
#    pedigree.add_individual('individual0', [1] * len(positions),genotype_likelihoods)
#    pedigree.add_individual('individual1', [1] * len(positions), [PhredGenotypeLikelihoods(1/3.0,1/3.0,1/3.0)] * len(positions))
#    pedigree.add_individual('individual2', [1] * len(positions), [PhredGenotypeLikelihoods(1/3.0,1/3.0,1/3.0)] * len(positions))
#    pedigree.add_relationship('individual0', 'individual1', 'individual2')
#    dp_forward_backward = GenotypeDPTable(numeric_sample_ids,readset,recombcost,pedigree)

# TODO
# check the results
#    compare_to_expected(dp_forward_backward, positions, expected, genotypes)

# first 5 tests compare to genotype likelihoods computed manually, with (non-)uniform priors


def test_geno_exact1():
    reads = """
          11
           01
        """

    expected_likelihoods = [
        PhredGenotypeLikelihoods([0.06666666666666667, 0.3333333333333333, 0.6]),
        PhredGenotypeLikelihoods([0.20930232558139536, 0.5813953488372093, 0.20930232558139536]),
        PhredGenotypeLikelihoods([0.06666666666666667, 0.3333333333333333, 0.6]),
    ]
    genotypes = canonic_index_list_to_biallelic_gt_list([2, 1, 2])
    check_genotyping_single_individual(reads, None, expected_likelihoods, genotypes, 10)


def test_geno_exact2():
    reads = """
        11
        11
        """
    weights = """
        11
        11
        """

    expected_likelihoods = [
        PhredGenotypeLikelihoods([0.00914139256727894, 0.25040580948312685, 0.7404527979495942]),
        PhredGenotypeLikelihoods([0.00914139256727894, 0.25040580948312685, 0.7404527979495942]),
    ]
    genotypes = canonic_index_list_to_biallelic_gt_list([2, 2])
    check_genotyping_single_individual(reads, weights, expected_likelihoods, genotypes, 10)


def test_geno_exact3():
    reads = """
          01
          11
        """

    expected_likelihoods = [
        PhredGenotypeLikelihoods([0.22163406214039125, 0.5567318757192175, 0.22163406214039125]),
        PhredGenotypeLikelihoods([0.009896432681242807, 0.18849252013808976, 0.8016110471806674]),
    ]
    check_genotyping_single_individual(reads, None, expected_likelihoods, None, 10)


def test_geno_priors1():
    reads = """
          01
          11
        """

    prior_likelihoods = [
        PhredGenotypeLikelihoods([0.1, 0.8, 0.1]),
        PhredGenotypeLikelihoods([0.1, 0.2, 0.7]),
    ]
    expected_likelihoods = [
        PhredGenotypeLikelihoods([0.04257892641700095, 0.9148421471659981, 0.04257892641700095]),
        PhredGenotypeLikelihoods([0.0016688611936185199, 0.05208684202468078, 0.9462442967817007]),
    ]
    check_genotyping_single_individual(
        reads, None, expected_likelihoods, None, 10, prior_likelihoods
    )


def test_geno_priors2():
    reads = """
            11
             01
             """

    prior_likelihoods = [
        PhredGenotypeLikelihoods([0, 0.5, 0.5]),
        PhredGenotypeLikelihoods([0.25, 0.5, 0.25]),
        PhredGenotypeLikelihoods([0.1, 0.4, 0.5]),
    ]
    expected_likelihoods = [
        PhredGenotypeLikelihoods([0.0, 0.35714285714285715, 0.6428571428571429]),
        PhredGenotypeLikelihoods([0.1323529411764706, 0.7352941176470589, 0.1323529411764706]),
        PhredGenotypeLikelihoods([0.015151515151515152, 0.30303030303030304, 0.6818181818181818]),
    ]
    check_genotyping_single_individual(
        reads, None, expected_likelihoods, None, 10, prior_likelihoods
    )


# check if exprected genotype predictions are made
def test_geno1():
    reads = """
    1111111111
    0000011111
    """

    genotypes = canonic_index_list_to_biallelic_gt_list([1, 1, 1, 1, 1, 2, 2, 2, 2, 2])
    check_genotyping_single_individual(reads, None, None, genotypes, 10)


def test_geno2():
    reads = """
    101
    101
    101
    101
    100
    100
    100
    100
    """

    genotypes = canonic_index_list_to_biallelic_gt_list([2, 0, 1])
    check_genotyping_single_individual(reads, None, None, genotypes, 10)


def test_geno3():
    reads = """
    111011
    110110
    110 10
    110110
    101110
    000 00
    01000
    000010
    100100
    """

    genotypes = canonic_index_list_to_biallelic_gt_list([1, 1, 0, 1, 1, 0])
    check_genotyping_single_individual(reads, None, None, genotypes, 10)


def test_geno4():
    reads = """
      1  11010
      00 00101
      001 01110
       1    111
    """

    check_genotyping_single_individual(reads, None, None, None, 10)


def test_geno5():
    reads = """
      0             0
      110111111111
      00100
           0001000000
           000
            10100
                  101
    """
    genotypes = canonic_index_list_to_biallelic_gt_list(
        [1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 0, 1]
    )
    check_genotyping_single_individual(reads, None, None, genotypes, 10)


def test_geno6():
    reads = """
        0100000000000
        0100010000000
        1110000000010
        0100000000000
        0101000001000
        0100010   000
        0 10000000100
        1111111011100
        0100111010011
        1111111000111
        1111110011111
        11110000  000
        1110011011111
        1111001011111
        0111111110  1
        """
    genotypes = canonic_index_list_to_biallelic_gt_list([1, 2, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1])
    check_genotyping_single_individual(reads, None, None, genotypes, 60)


def test_geno7():
    reads = """
        111
        101
        111
        101
        010
        000
        010
        000
        """
    genotypes = canonic_index_list_to_biallelic_gt_list([1, 1, 1])
    check_genotyping_single_individual(reads, None, None, genotypes, 60)


def test_geno8():
    reads = """
    11
    11
    10
    """
    genotypes = canonic_index_list_to_biallelic_gt_list([2, 1])
    check_genotyping_single_individual(reads, None, None, genotypes, 10)


def test_geno9():
    reads = """
    111
    111
    010
    010
       001
       001
       101
       101
    """
    genotypes = canonic_index_list_to_biallelic_gt_list([1, 2, 1, 1, 0, 2])
    check_genotyping_single_individual(reads, None, None, genotypes, 10)


def test_geno_10():
    reads = """
    001100
    000000
    000000
    110011
    110011
    111111
         """
    genotypes = canonic_index_list_to_biallelic_gt_list([1, 1, 0, 0, 1, 1])
    genotype_priors = [
        PhredGenotypeLikelihoods([0.1, 0.8, 0.1]),
        PhredGenotypeLikelihoods([0.1, 0.8, 0.1]),
        PhredGenotypeLikelihoods([0.7, 0.2, 0.1]),
        PhredGenotypeLikelihoods([0.7, 0.2, 0.1]),
        PhredGenotypeLikelihoods([0.1, 0.8, 0.1]),
        PhredGenotypeLikelihoods([0.1, 0.8, 0.1]),
    ]
    check_genotyping_single_individual(reads, None, None, genotypes, 50, genotype_priors)


def test_weighted_genotyping1():
    reads = """
      1  11010
      00 00101
      001 01110
       1    111
    """
    weights = """
      2  13112
      11 23359
      223 56789
       2    111
    """

    genotypes = canonic_index_list_to_biallelic_gt_list([1, 1, 2, 1, 1, 1, 2, 1, 1])
    check_genotyping_single_individual(reads, weights, None, genotypes, 10)


def test_weighted_genotyping2():
    reads = """
      111
      101
      111
      101
      010
      000
      010
      000
    """
    weights = """
      999
      999
      999
      999
      999
      999
      999
      999
    """

    # for the second position, each genotype should be equally likely
    expected_likelihoods = [
        PhredGenotypeLikelihoods([0, 1, 0]),
        PhredGenotypeLikelihoods([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]),
        PhredGenotypeLikelihoods([0, 1, 0]),
    ]
    check_genotyping_single_individual(reads, weights, expected_likelihoods, None, 50)


def test_weighted_genotyping3():
    reads = """
        0 1
         10
         """
    weights = """
        999
        999
    """
    expected_likelihoods = [
        PhredGenotypeLikelihoods([2.0 / 3.0, 1.0 / 3.0, 0]),
        PhredGenotypeLikelihoods([0, 1.0 / 3.0, 2.0 / 3.0]),
        PhredGenotypeLikelihoods([0, 1, 0]),
    ]
    check_genotyping_single_individual(reads, weights, expected_likelihoods, None, 300)


def test_weighted_genotyping4():
    reads = """
    00  00
    0000
    000
    111
    111101
    111111
    111110
      000
      1111
    """

    weights = """
    11  11
    1111
    111
    111
    111111
    111111
    111111
      111
      1111
    """
    genotypes = canonic_index_list_to_biallelic_gt_list([1, 1, 1, 1, 1, 1])
    check_genotyping_single_individual(reads, weights, None, genotypes, 10)


def test_weighted_genotyping5():
    reads = """
    1111
    1111
    1111
    1111
    1111
    1111
    1111
    1111
    1111
    1111
    1111
    0 00
    00
    0 00
    """

    weights = """
    1111
    1111
    1111
    1111
    1111
    1111
    1111
    1111
    1111
    1111
    1111
    1 11
    1111
    1 11
    """
    genotypes = canonic_index_list_to_biallelic_gt_list([1, 1, 1, 1])
    check_genotyping_single_individual(reads, weights, None, genotypes, 10)


def test_weighted_genotyping6():
    reads = """
        10
        10
         """
    weights = """
        99
        99
    """
    genotype_priors = [
        PhredGenotypeLikelihoods([0.5, 0.5, 0]),
        PhredGenotypeLikelihoods([0, 0.5, 0.5]),
    ]
    expected_likelihoods = [
        PhredGenotypeLikelihoods([0, 1, 0]),
        PhredGenotypeLikelihoods([0, 1, 0]),
    ]
    check_genotyping_single_individual(
        reads, weights, expected_likelihoods, None, 100, genotype_priors
    )


def test_small_example():
    reads = """
    11111111
    00000000
    """
    genotypes = canonic_index_list_to_biallelic_gt_list([1, 1, 1, 1, 1, 1, 1, 1])
    check_genotyping_single_individual(reads, None, None, genotypes, 1000)
