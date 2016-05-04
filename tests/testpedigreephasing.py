"""
Test phasing of pedigrees (PedMEC algorithm)
"""
from nose.tools import raises
from whatshap.core import PedigreeDPTable, ReadSet, Variant, Pedigree
from .phasingutils import string_to_readset, string_to_readset_pedigree, brute_force_phase


def trio_pedigree(*genotypes):
	"""
	order of genotypes: mother, father, child
	"""
	assert len(genotypes) == 3
	pedigree = Pedigree()
	pedigree.add_individual(0, genotypes[0])
	pedigree.add_individual(1, genotypes[1])
	pedigree.add_individual(2, genotypes[2])
	pedigree.add_relationship(0, 1, 2)
	return pedigree


def test_phase_empty_readset():
	rs = ReadSet()
	read_marks = []
	recombcost = []
	pedigree = trio_pedigree([], [], [])
	dp_table = PedigreeDPTable(rs, read_marks, recombcost, pedigree)
	(superreadsm, superreadsf, superreadsc), transmission_vector = dp_table.get_super_reads()


def compare_phasing_single_individual(reads, weights = None):
	"""Compares PedigreeDPTable based phasing to brute force phasing and returns string representation of superreads."""
	rs = string_to_readset(reads, weights)
	positions = rs.get_positions()
	read_marks = [0] * len(rs) # all reads from the child
	recombcost = [1] * len(positions) # recombination costs 1, should not occur 
	genotypesm = [1] * len(positions) # all genotypes heterozygous
	genotypesf = [1] * len(positions) # all genotypes heterozygous
	genotypesc = [1] * len(positions) # all genotypes heterozygous
	pedigree = trio_pedigree(genotypesm, genotypesf, genotypesc)
	dp_table = PedigreeDPTable(rs, read_marks, recombcost, pedigree)
	(superreadsm, superreadsf, superreadsc), transmission_vector = dp_table.get_super_reads()
	assert len(set(transmission_vector)) == 1
	assert len(superreadsc) == 2
	assert len(superreadsc[0]) == len(superreadsc[1])
	for v1, v2 in zip(*superreadsc):
		assert v1.position == v2.position
	haplotypes = tuple(sorted(''.join(str(v.allele) for v in sr) for sr in superreadsc))
	cost = dp_table.get_optimal_cost()
	partition = dp_table.get_optimal_partitioning()
	expected_cost, expected_partition, solution_count, expected_haplotype1, expected_haplotype2 = brute_force_phase(rs, True)
	inverse_partition = [1-p for p in partition]
	print()
	print(superreadsc[0])
	print(superreadsc[1])
	print('Partition:', partition)
	print('Expected: ', expected_partition)
	print('Haplotypes:')
	print(haplotypes[0])
	print(haplotypes[1])
	print('Expected Haplotypes:')
	print(expected_haplotype1)
	print(expected_haplotype2)
	print('Cost:', cost)
	print('Expected cost:', expected_cost)
	assert (partition == expected_partition) or (inverse_partition == expected_partition)
	assert solution_count == 1
	assert cost == expected_cost
	assert (haplotypes == (expected_haplotype1, expected_haplotype2)) or (haplotypes == (expected_haplotype2, expected_haplotype1))


def test_phase_single_individual_trivial() :
	reads = """
          11
           1
           01
        """
	compare_phasing_single_individual(reads)


def test_phase_single_individual1():
	reads = """
	 10
	 010
	 010
	"""


def test_phase_single_individual2():
	reads = """
	  1  11010
	  00 00101
	  001 0101
	"""
	compare_phasing_single_individual(reads)


def test_phase_single_individual3():
	reads = """
	  1  11010
	  00 00101
	  001 01010
	"""
	compare_phasing_single_individual(reads)


def test_phase_single_individual4():
	reads = """
	  1  11010
	  00 00101
	  001 01110
	   1    111
	"""
	compare_phasing_single_individual(reads)


def test_phase_single_individual4():
	reads = """
	  1  11010
	  00 00101
	  001 01110
	   1    111
	"""
	compare_phasing_single_individual(reads)


def test_phase_single_individual5():
	reads = """
	  0             0
	  110111111111
	  00100
	       0001000000
	       000
	        10100
	              101
	"""
	compare_phasing_single_individual(reads)


def test_weighted_phasing_single_individual1():
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
	compare_phasing_single_individual(reads, weights)


def phase_trio(reads, recombcost, genotypesm, genotypesf, genotypesc):
	rs, read_marks = string_to_readset_pedigree(reads)
	pedigree = trio_pedigree(genotypesm, genotypesf, genotypesc)
	dp_table = PedigreeDPTable(rs, read_marks, recombcost, pedigree)
	superreads_list, transmission_vector = dp_table.get_super_reads()
	cost = dp_table.get_optimal_cost()
	for superreads in superreads_list:
		for sr in superreads:
			print(sr)
	print('Cost:', dp_table.get_optimal_cost())
	print('Transmission vector:', transmission_vector)
	print('Partition:', dp_table.get_optimal_partitioning())
	return superreads_list, transmission_vector, cost


def assert_haplotypes(superreads_list, all_expected_haplotypes, length):
	for superreads, expected_haplotypes in zip(superreads_list, all_expected_haplotypes):
		assert len(superreads) == 2
		assert len(superreads[0]) == len(superreads[1]) == length
		haplotypes = tuple(sorted(''.join(str(v.allele) for v in sr) for sr in superreads))
		assert (haplotypes == (expected_haplotypes[0], expected_haplotypes[1])) or (haplotypes == (expected_haplotypes[1], expected_haplotypes[0]))


def test_phase_trio1() :
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
	genotypesm = [1,2,1]
	genotypesf = [1,1,1]
	genotypesc = [0,1,1]
	recombcost = [10,10,10]

	superreads_list, transmission_vector, cost = phase_trio(
		reads, recombcost, genotypesm, genotypesf, genotypesc)
	assert cost == 2
	assert len(set(transmission_vector)) == 1
	all_expected_haplotypes = [
		('111','010'),
		('001','110'),
		('001','010')
	]
	assert_haplotypes(superreads_list, all_expected_haplotypes, 3)


def test_phase_trio2() :
	reads = """
	  A 0
	  A 0
	  B 1
	  B 1
	  C 1
	  C 0
	"""
	genotypesm = [2]
	genotypesf = [0]
	genotypesc = [1]
	recombcost = [10,10,10]
	superreads_list, transmission_vector, cost = phase_trio(
		reads, recombcost, genotypesm, genotypesf, genotypesc)
	assert cost == 4
	assert len(set(transmission_vector)) == 1
	all_expected_haplotypes = [
		('1','1'),
		('0','0'),
		('0','1')
	]
	assert_haplotypes(superreads_list, all_expected_haplotypes, 1)


def test_phase_trio3() :
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
	genotypesm = [1,1,1,1,1,1]
	genotypesf = [1,1,1,1,1,1]
	genotypesc = [1,2,1,1,0,1]
	recombcost = [3,3,3,4,3,3]
	superreads_list, transmission_vector, cost = phase_trio(
		reads, recombcost, genotypesm, genotypesf, genotypesc)
	assert cost == 4
	assert transmission_vector in ([0,0,0,1,1,1], [1,1,1,0,0,0], [2,2,2,3,3,3], [3,3,3,2,2,2])
	all_expected_haplotypes = [
		('111111','000000'),
		('010101','101010'),
		('111000','010101')
	]
	assert_haplotypes(superreads_list, all_expected_haplotypes, 6)


def test_phase_trio4() :
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
	genotypesm = [1,1,1]
	genotypesf = [1,1,1]
	genotypesc = [1,1,1]
	recombcost = [1,1,1]
	superreads_list, transmission_vector, cost = phase_trio(
		reads, recombcost, genotypesm, genotypesf, genotypesc)
	assert cost == 2
	assert transmission_vector in ([0,2,0], [2,0,2], [1,3,1], [3,1,3])
	all_expected_haplotypes = [
		('111','000'),
		('101','010'),
		('111','000')
	]
	assert_haplotypes(superreads_list, all_expected_haplotypes, 3)


def test_phase_trio5() :
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
	genotypesm = [1,1,1]
	genotypesf = [1,1,1]
	genotypesc = [1,1,1]
	recombcost = [2,2,2]
	superreads_list, transmission_vector, cost = phase_trio(
		reads, recombcost, genotypesm, genotypesf, genotypesc)
	assert cost == 3
	assert len(set(transmission_vector)) == 1
	all_expected_haplotypes = [
		('111','000'),
		('111','000'),
		('111','000')
	]
	assert_haplotypes(superreads_list, all_expected_haplotypes, 3)
