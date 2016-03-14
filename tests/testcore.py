from nose.tools import raises
from whatshap.core import Read, DPTable, ReadSet, Variant
from .phasingutils import string_to_readset, string_to_readset_trio, brute_force_phase


def test_read():
	r = Read("name", 15)
	assert r.name == "name"
	assert r.mapqs[0] == 15

	assert r.is_sorted

	r.add_variant(100, 1, 37)
	r.add_variant(23, 0, 99)
	assert not r.is_sorted
	r.sort()
	assert r.is_sorted

	assert 100 in r
	assert 23 in r
	assert not 22 in r
	assert not 24 in r
	assert not 1000 in r
	assert not -1000 in r


def test_read_iteration():
	r = Read("name", 15)
	r.add_variant(100, 1, 37)
	r.add_variant(23, 0, 99)
	v1 = Variant(position=100, allele=1, quality=37)
	v2 = Variant(position=23, allele=0, quality=99)
	variants = list(r)
	assert variants == [v1, v2]
	# negative indices
	assert r[-1] == v2
	assert r[-2] == v1


@raises(IndexError)
def test_read_indexerror1():
	r = Read("name", 15)
	r.add_variant(100, 1, 37)
	r.add_variant(23, 0, 99)
	r[2]


@raises(IndexError)
def test_read_indexerror2():
	r = Read("name", 15)
	r.add_variant(100, 1, 37)
	r.add_variant(23, 0, 99)
	r[-3]


def test_empty_readset():
	rs = ReadSet()
	assert len(rs) == 0


def test_readset():
	rs = ReadSet()
	r = Read('Read A', 56)
	r.add_variant(100, 1, 37)
	r.add_variant(101, 0, 18)
	rs.add(r)

	r = Read('Read B', 0)
	r.add_variant(101, 0, 23)
	rs.add(r)

	r = Read('Read C', 17)
	r.add_variant(99, 1, 27)
	r.add_variant(80, 1, 17)
	r[1] = Variant(position=105, allele=0, quality=14)
	rs.add(r)

	assert rs[0].name == 'Read A'
	assert rs[1].name == 'Read B'
	assert rs[2].name == 'Read C'

	rs.sort()

	# should be sorted after finalization
	assert rs[0].name == 'Read C'
	assert rs[1].name == 'Read A'
	assert rs[2].name == 'Read B'

	assert len(rs) == 3

	assert rs.get_positions() == [99, 100, 101, 105]

	r = rs[(0,'Read A')]
	assert r.name == 'Read A'
	assert r.mapqs == (56,), str(r.mapqs)

	r = rs[(0,'Read B')]
	assert r.name == 'Read B'
	assert r.mapqs == (0,)

	r = rs[(0,'Read C')]
	assert r.name == 'Read C'
	assert r.mapqs == (17,)
	assert len(r) == 2
	assert r[0] == Variant(position=99, allele=1, quality=27)
	assert r[1] == Variant(position=105, allele=0, quality=14)

def test_readset2():
	rs = ReadSet()
	rs.add(Read('Read A', 1, 23))
	rs.add(Read('Read A', 2, 70))
	rs.add(Read('Read B', 3, 23))
	assert rs[(23,'Read A')].mapqs == (1,)
	assert rs[(70,'Read A')].mapqs == (2,)
	assert rs[(23,'Read B')].mapqs == (3,)


@raises(KeyError)
def test_non_existing_read_name():
	rs = ReadSet()
	r = Read('Read A', 56)
	r.add_variant(100, 1, 37)
	r.add_variant(101, 0, 18)
	rs.add(r)
	rs[(0,'foo')]


@raises(KeyError)
def test_non_existing_read_name2():
	rs = ReadSet()
	r = Read('Read A', 56, 1)
	r.add_variant(100, 1, 37)
	r.add_variant(101, 0, 18)
	rs.add(r)
	rs[(2,'Read A')]


# TODO: Test subset method


def test_phase_empty_readset():
	rs = ReadSet()
	read_marks = []
	recombcost = []
	genotypesm, genotypesf, genotypesc = [], [], []
	dp_table = DPTable(rs, read_marks, recombcost, genotypesm, genotypesf, genotypesc)
	superreadsm, superreadsf, superreadsc, transmission_vector = dp_table.get_super_reads()


def compare_phasing_single_individual(reads, weights = None):
	"""Compares DPTable based phasing to brute force phasing and returns string representation of superreads."""
	rs = string_to_readset(reads, weights)
	positions = rs.get_positions()
	read_marks = [0] * len(rs) # all reads from the child
	recombcost = [1] * len(positions) # recombination costs 1, should not occur 
	genotypesm = [1] * len(positions) # all genotypes heterozygous
	genotypesf = [1] * len(positions) # all genotypes heterozygous
	genotypesc = [1] * len(positions) # all genotypes heterozygous
	dp_table = DPTable(rs, read_marks, recombcost, genotypesm, genotypesf, genotypesc)
	superreadsm, superreadsf, superreadsc, transmission_vector = dp_table.get_super_reads()
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


def test_phase_trio1() :
	reads = """
	  M 111
	  M 010
	  M 110
	  F 001
	  F 110
	  F 101
	  C 001
	  C 010
	  C 010
	"""
	genotypesm = [1,2,1]
	genotypesf = [1,1,1]
	genotypesc = [0,1,1]
	recombcost = [10,10,10]
	rs, read_marks = string_to_readset_trio(reads)
	dp_table = DPTable(rs, read_marks, recombcost, genotypesm, genotypesf, genotypesc)
	assert dp_table.get_optimal_cost() == 2
	superreadsm, superreadsf, superreadsc, transmission_vector = dp_table.get_super_reads()
	assert len(set(transmission_vector)) == 1
	all_expected_haplotypes = [
		('111','010'),
		('001','110'),
		('001','010')
	]
	for superreads, expected_haplotypes in zip([superreadsm, superreadsf, superreadsc],all_expected_haplotypes):
		assert len(superreads) == 2
		assert len(superreads[0]) == len(superreads[1]) == 3
		haplotypes = tuple(sorted(''.join(str(v.allele) for v in sr) for sr in superreads))
		assert (haplotypes == (expected_haplotypes[0], expected_haplotypes[1])) or (haplotypes == (expected_haplotypes[1], expected_haplotypes[0]))


def test_phase_trio2() :
	reads = """
	  M 0
	  M 0
	  F 1
	  F 1
	  C 1
	  C 0
	"""
	genotypesm = [2]
	genotypesf = [0]
	genotypesc = [1]
	recombcost = [10,10,10]
	rs, read_marks = string_to_readset_trio(reads)
	dp_table = DPTable(rs, read_marks, recombcost, genotypesm, genotypesf, genotypesc)
	assert dp_table.get_optimal_cost() == 4
	superreadsm, superreadsf, superreadsc, transmission_vector = dp_table.get_super_reads()
	assert len(set(transmission_vector)) == 1
	all_expected_haplotypes = [
		('1','1'),
		('0','0'),
		('0','1')
	]
	for superreads in [superreadsm, superreadsf, superreadsc]:
		for sr in superreads:
			print(sr)
	for superreads, expected_haplotypes in zip([superreadsm, superreadsf, superreadsc],all_expected_haplotypes):
		assert len(superreads) == 2
		assert len(superreads[0]) == len(superreads[1]) == 1
		haplotypes = tuple(sorted(''.join(str(v.allele) for v in sr) for sr in superreads))
		assert (haplotypes == (expected_haplotypes[0], expected_haplotypes[1])) or (haplotypes == (expected_haplotypes[1], expected_haplotypes[0]))


def test_phase_trio3() :
	reads = """
	  M 1111
	  F 1010
	  C 111000
	  C 010101
	  F 0101
	  M  0000
	  F  1010
	  C  1010
	  C  1100
	  M   0000
	  M   1111
	  F   1010
	  F    010
	"""
	genotypesm = [1,1,1,1,1,1]
	genotypesf = [1,1,1,1,1,1]
	genotypesc = [1,2,1,1,0,1]
	recombcost = [4,4,4,4,4,4]
	rs, read_marks = string_to_readset_trio(reads)
	dp_table = DPTable(rs, read_marks, recombcost, genotypesm, genotypesf, genotypesc)
	superreadsm, superreadsf, superreadsc, transmission_vector = dp_table.get_super_reads()
	for superreads in [superreadsm, superreadsf, superreadsc]:
		for sr in superreads:
			print(sr)
	print('Cost:',dp_table.get_optimal_cost())
	print('Transmission vector:', transmission_vector)
	assert dp_table.get_optimal_cost() == 4
	assert transmission_vector in ([0,0,0,1,1,1], [1,1,1,0,0,0], [2,2,2,3,3,3], [3,3,3,2,2,2])
	all_expected_haplotypes = [
		('111111','000000'),
		('010101','101010'),
		('111000','010101')
	]
	for superreads, expected_haplotypes in zip([superreadsm, superreadsf, superreadsc],all_expected_haplotypes):
		assert len(superreads) == 2
		assert len(superreads[0]) == len(superreads[1]) == 6
		haplotypes = tuple(sorted(''.join(str(v.allele) for v in sr) for sr in superreads))
		assert (haplotypes == (expected_haplotypes[0], expected_haplotypes[1])) or (haplotypes == (expected_haplotypes[1], expected_haplotypes[0]))
