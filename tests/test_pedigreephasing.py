"""
Test phasing of pedigrees (PedMEC algorithm)
"""
from collections import defaultdict
from pytest import raises
from whatshap.core import (
    PedigreeDPTable,
    ReadSet,
    Pedigree,
    NumericSampleIds,
    PhredGenotypeLikelihoods,
)
from whatshap.pedigree import centimorgen_to_phred
from whatshap.testhelpers import string_to_readset_pedigree, canonic_index_list_to_biallelic_gt_list


def phase_pedigree(reads, recombcost, pedigree, distrust_genotypes=False, positions=None):
    rs = string_to_readset_pedigree(reads)
    dp_table = PedigreeDPTable(rs, recombcost, pedigree, distrust_genotypes, positions)
    superreads_list, transmission_vector = dp_table.get_super_reads()
    cost = dp_table.get_optimal_cost()
    for superreads in superreads_list:
        for sr in superreads:
            print(sr)
    print("Cost:", dp_table.get_optimal_cost())
    print("Transmission vector:", transmission_vector)
    print("Partition:", dp_table.get_optimal_partitioning())
    return superreads_list, transmission_vector, cost


def assert_haplotypes(superreads_list, all_expected_haplotypes, length):
    for superreads, expected_haplotypes in zip(superreads_list, all_expected_haplotypes):
        assert len(superreads) == 2
        assert len(superreads[0]) == len(superreads[1]) == length
        haplotypes = tuple(sorted("".join(str(v.allele) for v in sr) for sr in superreads))
        assert (haplotypes == (expected_haplotypes[0], expected_haplotypes[1])) or (
            haplotypes == (expected_haplotypes[1], expected_haplotypes[0])
        )


def assert_trio_allele_order(superreads_list, transmission_vector, nr_of_positions):
    # assume superreads_list contains superreads for father, mother, child (in that order!)
    assert len(superreads_list) == 3
    father = superreads_list[0]
    mother = superreads_list[1]
    child = superreads_list[2]

    for pos in range(nr_of_positions):
        transmission_value = transmission_vector[pos]
        paternal_transmission = transmission_value % 2
        maternal_transmission = transmission_value // 2
        paternal_allele = father[not paternal_transmission][pos].allele
        maternal_allele = mother[not maternal_transmission][pos].allele
        child_allele_p = child[0][pos].allele
        child_allele_m = child[1][pos].allele
        print(
            "position: ",
            pos,
            "paternal allele: ",
            paternal_allele,
            "maternal allele",
            maternal_allele,
            "child genotype: ",
            child_allele_p,
            child_allele_m,
        )
        assert paternal_allele == child_allele_p
        assert maternal_allele == child_allele_m


def get_trio_transmission_vectors(transmission_vector, nr_of_trios):
    trio_transmission_vectors = defaultdict(list)
    for transmission_value in transmission_vector:
        for trio in range(nr_of_trios):
            value = transmission_value % 4
            transmission_value = transmission_value // 4
            trio_transmission_vectors[trio].append(value)
    return trio_transmission_vectors


def test_phase_empty_trio():
    rs = ReadSet()
    recombcost = []
    pedigree = Pedigree(NumericSampleIds())
    pedigree.add_individual("individual0", [])
    pedigree.add_individual("individual1", [])
    pedigree.add_individual("individual2", [])
    pedigree.add_relationship("individual0", "individual1", "individual2")
    dp_table = PedigreeDPTable(rs, recombcost, pedigree)
    ((superreadsm, superreadsf, superreadsc), transmission_vector) = dp_table.get_super_reads()


def test_phase_trio1():
    reads = """
      A 111
      A 010
      A 110
      B 001
      B 110
      B 101
      C 001
      C 010
      C 010
    """
    pedigree = Pedigree(NumericSampleIds())
    pedigree.add_individual("individual0", canonic_index_list_to_biallelic_gt_list([1, 2, 1]))
    pedigree.add_individual("individual1", canonic_index_list_to_biallelic_gt_list([1, 1, 1]))
    pedigree.add_individual("individual2", canonic_index_list_to_biallelic_gt_list([0, 1, 1]))
    pedigree.add_relationship("individual0", "individual1", "individual2")
    recombcost = [10, 10, 10]
    superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
    assert cost == 2
    assert len(set(transmission_vector)) == 1
    all_expected_haplotypes = [("111", "010"), ("001", "110"), ("010", "001")]
    assert_haplotypes(superreads_list, all_expected_haplotypes, 3)
    assert_trio_allele_order(superreads_list, transmission_vector, 3)


def test_phase_trio2():
    reads = """
      A 00
      A 00
      B 11
      B 11
      C 11
      C 00
    """
    pedigree = Pedigree(NumericSampleIds())
    pedigree.add_individual("individual0", canonic_index_list_to_biallelic_gt_list([2, 2]))
    pedigree.add_individual("individual1", canonic_index_list_to_biallelic_gt_list([0, 0]))
    pedigree.add_individual("individual2", canonic_index_list_to_biallelic_gt_list([1, 1]))
    pedigree.add_relationship("individual0", "individual1", "individual2")
    recombcost = [10, 10, 10]
    superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
    assert cost == 8
    assert len(set(transmission_vector)) == 1
    all_expected_haplotypes = [("11", "11"), ("00", "00"), ("00", "11")]
    assert_haplotypes(superreads_list, all_expected_haplotypes, 2)
    assert_trio_allele_order(superreads_list, transmission_vector, 2)


def test_phase_trio3():
    reads = """
      A 1111
      B 1010
      C 111000
      C 010101
      B 0101
      A  0000
      B  1010
      C  1010
      C  1100
      A   0000
      A   1111
      B   1010
      B    010
    """
    pedigree = Pedigree(NumericSampleIds())
    pedigree.add_individual(
        "individual0", canonic_index_list_to_biallelic_gt_list([1, 1, 1, 1, 1, 1])
    )
    pedigree.add_individual(
        "individual1", canonic_index_list_to_biallelic_gt_list([1, 1, 1, 1, 1, 1])
    )
    pedigree.add_individual(
        "individual2", canonic_index_list_to_biallelic_gt_list([1, 2, 1, 1, 0, 1])
    )
    pedigree.add_relationship("individual0", "individual1", "individual2")
    recombcost = [3, 3, 3, 4, 3, 3]
    superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
    assert cost == 4
    assert transmission_vector in (
        [0, 0, 0, 1, 1, 1],
        [1, 1, 1, 0, 0, 0],
        [2, 2, 2, 3, 3, 3],
        [3, 3, 3, 2, 2, 2],
    )
    all_expected_haplotypes = [
        ("111111", "000000"),
        ("010101", "101010"),
        ("111000", "010101"),
    ]
    assert_haplotypes(superreads_list, all_expected_haplotypes, 6)
    assert_trio_allele_order(superreads_list, transmission_vector, 6)


def test_phase_trio4():
    reads = """
      B 101
      B 101
      B 101
      A 111
      A 111
      A 111
      C 111
      C 111
      C 111
    """
    pedigree = Pedigree(NumericSampleIds())
    pedigree.add_individual("individual0", canonic_index_list_to_biallelic_gt_list([1, 1, 1]))
    pedigree.add_individual("individual1", canonic_index_list_to_biallelic_gt_list([1, 1, 1]))
    pedigree.add_individual("individual2", canonic_index_list_to_biallelic_gt_list([1, 1, 1]))
    pedigree.add_relationship("individual0", "individual1", "individual2")
    recombcost = [1, 1, 1]
    superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
    assert cost == 2
    assert transmission_vector in ([0, 2, 0], [2, 0, 2], [1, 3, 1], [3, 1, 3])
    all_expected_haplotypes = [("111", "000"), ("101", "010"), ("111", "000")]
    assert_haplotypes(superreads_list, all_expected_haplotypes, 3)
    assert_trio_allele_order(superreads_list, transmission_vector, 3)


def test_phase_trio5():
    reads = """
      B 101
      B 101
      B 101
      A 111
      A 111
      A 111
      C 111
      C 111
      C 111
    """
    pedigree = Pedigree(NumericSampleIds())
    pedigree.add_individual("individual0", canonic_index_list_to_biallelic_gt_list([1, 1, 1]))
    pedigree.add_individual("individual1", canonic_index_list_to_biallelic_gt_list([1, 1, 1]))
    pedigree.add_individual("individual2", canonic_index_list_to_biallelic_gt_list([1, 1, 1]))
    pedigree.add_relationship("individual0", "individual1", "individual2")
    recombcost = [2, 2, 2]
    superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
    assert cost == 3
    assert len(set(transmission_vector)) == 1
    all_expected_haplotypes = [("111", "000"), ("111", "000"), ("111", "000")]
    assert_haplotypes(superreads_list, all_expected_haplotypes, 3)
    assert_trio_allele_order(superreads_list, transmission_vector, 3)


def test_phase_trio_pure_genetic():
    reads = ""
    pedigree = Pedigree(NumericSampleIds())
    pedigree.add_individual("individual0", canonic_index_list_to_biallelic_gt_list([2, 1, 1, 0]))
    pedigree.add_individual("individual1", canonic_index_list_to_biallelic_gt_list([1, 2, 2, 1]))
    pedigree.add_individual("individual2", canonic_index_list_to_biallelic_gt_list([1, 1, 1, 0]))
    pedigree.add_relationship("individual0", "individual1", "individual2")
    recombcost = [2, 2, 2]
    superreads_list, transmission_vector, cost = phase_pedigree(
        reads, recombcost, pedigree, positions=[10, 20, 30, 40]
    )
    assert cost == 0
    assert len(set(transmission_vector)) == 1
    all_expected_haplotypes = [("1110", "1000"), ("1111", "0110"), ("1000", "0110")]
    assert_haplotypes(superreads_list, all_expected_haplotypes, 4)
    assert_trio_allele_order(superreads_list, transmission_vector, 4)


def test_phase_doubletrio_pure_genetic():
    reads = ""
    pedigree = Pedigree(NumericSampleIds())
    pedigree.add_individual("individualA", canonic_index_list_to_biallelic_gt_list([1, 2, 1, 0]))
    pedigree.add_individual("individualB", canonic_index_list_to_biallelic_gt_list([1, 0, 1, 1]))
    pedigree.add_individual("individualC", canonic_index_list_to_biallelic_gt_list([2, 1, 1, 0]))
    pedigree.add_individual("individualD", canonic_index_list_to_biallelic_gt_list([1, 2, 2, 1]))
    pedigree.add_individual("individualE", canonic_index_list_to_biallelic_gt_list([1, 1, 1, 0]))
    pedigree.add_relationship("individualA", "individualB", "individualC")
    pedigree.add_relationship("individualC", "individualD", "individualE")
    recombcost = [2, 2, 2]
    superreads_list, transmission_vector, cost = phase_pedigree(
        reads, recombcost, pedigree, positions=[10, 20, 30, 40]
    )
    assert cost == 0
    assert len(set(transmission_vector)) == 1
    all_expected_haplotypes = [
        ("0100", "1110"),
        ("0011", "1000"),
        ("1110", "1000"),
        ("1111", "0110"),
        ("1000", "0110"),
    ]
    assert_haplotypes(superreads_list, all_expected_haplotypes, 4)
    trio_transmission_vectors = get_trio_transmission_vectors(transmission_vector, 4)
    assert_trio_allele_order(superreads_list[:3], trio_transmission_vectors[0], 4)
    assert_trio_allele_order(superreads_list[2:], trio_transmission_vectors[1], 4)


def test_phase_quartet1():
    reads = """
      A 111
      A 010
      A 110
      B 001
      B 110
      B 101
      C 001
      C 010
      C 010
      D 001
      D 010
      D 010
    """
    pedigree = Pedigree(NumericSampleIds())
    pedigree.add_individual("individual0", canonic_index_list_to_biallelic_gt_list([1, 2, 1]))
    pedigree.add_individual("individual1", canonic_index_list_to_biallelic_gt_list([1, 1, 1]))
    pedigree.add_individual("individual2", canonic_index_list_to_biallelic_gt_list([0, 1, 1]))
    pedigree.add_individual("individual3", canonic_index_list_to_biallelic_gt_list([0, 1, 1]))
    pedigree.add_relationship("individual0", "individual1", "individual2")
    pedigree.add_relationship("individual0", "individual1", "individual3")
    recombcost = [10, 10, 10]
    superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
    assert cost == 2
    assert len(set(transmission_vector)) == 1
    all_expected_haplotypes = [
        ("111", "010"),
        ("001", "110"),
        ("001", "010"),
        ("001", "010"),
    ]
    assert_haplotypes(superreads_list, all_expected_haplotypes, 3)
    trio_transmission_vectors = get_trio_transmission_vectors(transmission_vector, 3)
    assert_trio_allele_order(superreads_list[:3], trio_transmission_vectors[0], 3)
    assert_trio_allele_order(
        [superreads_list[0], superreads_list[1], superreads_list[3]],
        trio_transmission_vectors[1],
        3,
    )


def test_phase_quartet2():
    reads = """
      A 111111
      A 000000
      B 010101
      B 101010
      C 000000
      C 010101
      D 000000
      D 010101
    """
    pedigree = Pedigree(NumericSampleIds())
    pedigree.add_individual(
        "individual0", canonic_index_list_to_biallelic_gt_list([1, 1, 1, 1, 1, 1])
    )
    pedigree.add_individual(
        "individual1", canonic_index_list_to_biallelic_gt_list([1, 1, 1, 1, 1, 1])
    )
    pedigree.add_individual(
        "individual2", canonic_index_list_to_biallelic_gt_list([0, 1, 0, 1, 0, 1])
    )
    pedigree.add_individual(
        "individual3", canonic_index_list_to_biallelic_gt_list([0, 1, 0, 1, 0, 1])
    )
    pedigree.add_relationship("individual0", "individual1", "individual2")
    pedigree.add_relationship("individual0", "individual1", "individual3")
    recombcost = [3, 3, 3, 3, 3, 3]

    superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
    assert cost == 0
    assert len(set(transmission_vector)) == 1
    all_expected_haplotypes = [
        ("111111", "000000"),
        ("010101", "101010"),
        ("000000", "010101"),
        ("000000", "010101"),
    ]
    assert_haplotypes(superreads_list, all_expected_haplotypes, 6)
    trio_transmission_vectors = get_trio_transmission_vectors(transmission_vector, 6)
    assert_trio_allele_order(superreads_list[:3], trio_transmission_vectors[0], 6)
    assert_trio_allele_order(
        [superreads_list[0], superreads_list[1], superreads_list[3]],
        trio_transmission_vectors[1],
        6,
    )


def test_phase_quartet3():
    reads = """
      A 1111
      A 0000
      B 1010
      C 111000
      C 010101
      D 000000
      D 010
      B 0101
      C  1100
      D  10010
      A   0000
      A   1111
      B   1010
      B   0101
    """
    pedigree = Pedigree(NumericSampleIds())
    pedigree.add_individual(
        "individual0", canonic_index_list_to_biallelic_gt_list([1, 1, 1, 1, 1, 1])
    )
    pedigree.add_individual(
        "individual1", canonic_index_list_to_biallelic_gt_list([1, 1, 1, 1, 1, 1])
    )
    pedigree.add_individual(
        "individual2", canonic_index_list_to_biallelic_gt_list([1, 2, 1, 1, 0, 1])
    )
    pedigree.add_individual(
        "individual3", canonic_index_list_to_biallelic_gt_list([0, 1, 0, 0, 1, 0])
    )
    pedigree.add_relationship("individual0", "individual1", "individual2")
    pedigree.add_relationship("individual0", "individual1", "individual3")
    recombcost = [3, 3, 3, 4, 3, 3]
    superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree)
    print(cost)
    print(transmission_vector)
    assert cost == 8
    # TODO: expect transmission in both trio relations. Update once transmission vectors
    #       are returned per trio relationship
    # assert transmission_vector in ([0,0,0,1,1,1], [1,1,1,0,0,0], [2,2,2,3,3,3], [3,3,3,2,2,2])
    all_expected_haplotypes = [
        ("111111", "000000"),
        ("010101", "101010"),
        ("111000", "010101"),
        ("000000", "010010"),
    ]
    assert_haplotypes(superreads_list, all_expected_haplotypes, 6)
    trio_transmission_vectors = get_trio_transmission_vectors(transmission_vector, 6)
    assert_trio_allele_order(superreads_list[:3], trio_transmission_vectors[0], 6)
    assert_trio_allele_order(
        [superreads_list[0], superreads_list[1], superreads_list[3]],
        trio_transmission_vectors[1],
        6,
    )


def test_centimorgen_to_phred():
    assert round(centimorgen_to_phred(0.10010013353365396)) == 30
    assert round(centimorgen_to_phred(0.0010000100001343354)) == 50
    assert round(centimorgen_to_phred(1e-38)) == 400


def test_centimorgen_to_phred_zero():
    with raises(ValueError):
        assert centimorgen_to_phred(0)


def test_phase_trio_genotype_likelihoods():
    reads = """
      A 111
      A 010
      A 110
      B 001
      B 110
      B 101
      C 001
      C 010
      C 010
    """
    pedigree = Pedigree(NumericSampleIds())
    genotype_likelihoods_mother = [
        PhredGenotypeLikelihoods([0, 0, 0]),
        PhredGenotypeLikelihoods([0, 0, 1]),
        PhredGenotypeLikelihoods([5, 0, 5]),
    ]
    genotype_likelihoods0 = [PhredGenotypeLikelihoods([0, 0, 0])] * 3
    pedigree.add_individual(
        "individual0",
        canonic_index_list_to_biallelic_gt_list([0, 0, 0]),
        genotype_likelihoods_mother,
    )
    pedigree.add_individual(
        "individual1", canonic_index_list_to_biallelic_gt_list([0, 0, 0]), genotype_likelihoods0
    )
    pedigree.add_individual(
        "individual2", canonic_index_list_to_biallelic_gt_list([0, 0, 0]), genotype_likelihoods0
    )
    pedigree.add_relationship("individual0", "individual1", "individual2")
    recombcost = [10, 10, 10]
    superreads_list, transmission_vector, cost = phase_pedigree(reads, recombcost, pedigree, True)
    assert cost == 3
    assert len(set(transmission_vector)) == 1
    all_expected_haplotypes = [("111", "010"), ("001", "110"), ("001", "010")]
    assert_haplotypes(superreads_list, all_expected_haplotypes, 3)
    assert_trio_allele_order(superreads_list, transmission_vector, 3)
