from nose.tools import raises
from whatshap.core import DPTable, ReadSet
from .phasingutils import string_to_readset, brute_force_phase


def test_phase_empty_readset():
	rs = ReadSet()
	dp_table = DPTable(rs, all_heterozygous=False)
	superreads = dp_table.get_super_reads()


def compare_phasing(reads, all_heterozygous, weights = None):
	"""Compares DPTable based phasing to brute force phasing and returns string representation of superreads."""
	rs = string_to_readset(reads, weights)
	dp_table = DPTable(rs, all_heterozygous)
	superreads = dp_table.get_super_reads()
	assert len(superreads) == 2
	assert len(superreads[0]) == len(superreads[1])
	for v1, v2 in zip(*superreads):
		assert v1.position == v2.position
	haplotypes = tuple(sorted(''.join(str(v.allele) for v in sr) for sr in superreads))
	cost = dp_table.get_optimal_cost()
	partition = dp_table.get_optimal_partitioning()
	expected_cost, expected_partition, solution_count, expected_haplotype1, expected_haplotype2 = brute_force_phase(rs, all_heterozygous)
	inverse_partition = [1-p for p in partition]
	print()
	print(superreads[0])
	print(superreads[1])
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


def test_phase_trivial() :
	reads = """
          11
           1
           01
        """
	compare_phasing(reads, True)
	compare_phasing(reads, False)


def test_phase1():
	reads = """
	 10
	 010
	 010
	"""
	compare_phasing(reads, True)
	compare_phasing(reads, False)


def test_phase2():
	reads = """
	  1  11010
	  00 00101
	  001 0101
	"""
	compare_phasing(reads, True)
	compare_phasing(reads, False)


def test_phase3():
	reads = """
	  1  11010
	  00 00101
	  001 01010
	"""
	compare_phasing(reads, True)
	compare_phasing(reads, False)


def test_phase4():
	reads = """
	  1  11010
	  00 00101
	  001 01110
	   1    111
	"""
	compare_phasing(reads, True)
	compare_phasing(reads, False)


def test_phase4():
	reads = """
	  1  11010
	  00 00101
	  001 01110
	   1    111
	"""
	compare_phasing(reads, True)
	compare_phasing(reads, False)


def test_phase5():
	reads = """
	  0             0
	  110111111111
	  00100
	       0001000000
	       000
	        10100
	              101
	"""
	compare_phasing(reads, True)
	compare_phasing(reads, False)


def test_weighted_phasing1():
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
	compare_phasing(reads, True, weights)
	compare_phasing(reads, False, weights)
