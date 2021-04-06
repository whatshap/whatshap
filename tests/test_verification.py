from whatshap.core import PedigreeDPTable, Pedigree, NumericSampleIds, PhredGenotypeLikelihoods
from whatshap.testhelpers import string_to_readset, matrix_to_readset, canonic_index_to_biallelic_gt
from whatshap.verification import verify_mec_score_and_partitioning


def verify(rs, all_heterozygous=False):
    positions = rs.get_positions()
    # recombination costs 1, should not occur
    recombcost = [1] * len(positions)
    pedigree = Pedigree(NumericSampleIds())
    genotype_likelihoods = [
        None if all_heterozygous else PhredGenotypeLikelihoods([0, 0, 0])
    ] * len(positions)
    # all genotypes heterozygous
    pedigree.add_individual(
        "individual0",
        [canonic_index_to_biallelic_gt(1) for _ in range(len(positions))],
        genotype_likelihoods,
    )
    dp_table = PedigreeDPTable(rs, recombcost, pedigree, distrust_genotypes=not all_heterozygous)
    verify_mec_score_and_partitioning(dp_table, rs)


def test_string():
    reads = """
      0             0
      110111111111
      00100
           0001000000
           000
            10100
                  101
    """
    rs = string_to_readset(reads)
    verify(rs, True)
    verify(rs, False)


def test_matrix():
    with open("tests/test.matrix") as f:
        rs = matrix_to_readset(f)
    verify(rs, True)
    verify(rs, False)
