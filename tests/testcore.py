from nose.tools import raises
from whatshap.core import Read, DPTable, ReadSet, Variant
from .phasingutils import string_to_readset, brute_force_phase


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
