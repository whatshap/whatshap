import textwrap
from nose.tools import raises
from whatshap.core import Read, DPTable, ReadSet, Variant


def test_read():
	r = Read("name", 15)
	assert r.name == "name"
	assert r.mapqs[0] == 15

	assert r.is_sorted

	r.add_variant(100, 'A', 1, 37)
	r.add_variant(23, 'T', 0, 99)
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
	r.add_variant(100, 'A', 1, 37)
	r.add_variant(23, 'T', 0, 99)
	v1 = Variant(position=100, base='A', allele=1, quality=37)
	v2 = Variant(position=23, base='T', allele=0, quality=99)
	variants = list(r)
	assert variants == [v1, v2]
	# negative indices
	assert r[-1] == v2
	assert r[-2] == v1


@raises(IndexError)
def test_read_indexerror1():
	r = Read("name", 15)
	r.add_variant(100, 'A', 1, 37)
	r.add_variant(23, 'T', 0, 99)
	r[2]


@raises(IndexError)
def test_read_indexerror2():
	r = Read("name", 15)
	r.add_variant(100, 'A', 1, 37)
	r.add_variant(23, 'T', 0, 99)
	r[-3]


def test_empty_readset():
	rs = ReadSet()
	assert len(rs) == 0


def test_readset():
	rs = ReadSet()
	r = Read('Read A', 56)
	r.add_variant(100, 'A', 1, 37)
	r.add_variant(101, 'C', 0, 18)
	rs.add(r)

	r = Read('Read B', 0)
	r.add_variant(101, 'C', 0, 23)
	rs.add(r)

	r = Read('Read C', 17)
	r.add_variant(99, 'G', 1, 27)
	r.add_variant(105, 'T', 0, 14)
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

	r = rs['Read A']
	assert r.name == 'Read A'
	assert r.mapqs == (56,), str(r.mapqs)

	r = rs['Read B']
	assert r.name == 'Read B'
	assert r.mapqs == (0,)


@raises(KeyError)
def test_non_existing_read_name():
	rs = ReadSet()
	r = Read('Read A', 56)
	r.add_variant(100, 'A', 1, 37)
	r.add_variant(101, 'C', 0, 18)
	rs.add(r)
	rs['foo']


# TODO: Test subset method


def test_phase_empty_readset():
	rs = ReadSet()
	dp_table = DPTable(rs, all_heterozygous=False)
	superreads = dp_table.get_super_reads()


def string_to_readset(s):
	s = textwrap.dedent(s).strip()
	bits = { '0': 'A', '1': 'C', 'E': 'G' }
	rs = ReadSet()
	for index, line in enumerate(s.split('\n'), 1):
		read = Read('Read {}'.format(index), 50)
		for pos, c in enumerate(line, 1):
			if c == ' ':
				continue
			read.add_variant(position=pos * 10, base=bits[c], allele=int(c), quality=1)
		rs.add(read)
	print(rs)
	return rs


def flip_cost(variant, target_value):
	"""Returns cost of flipping the given read variant to target_value."""
	if variant.allele == target_value:
		return 0
	else:
		return variant.quality


def column_cost(variants, possible_assignments):
	"""Compute cost for one position."""
	costs = []
	for allele1, allele2 in possible_assignments:
		cost1 = sum(flip_cost(v,allele1) for v in variants[0])
		cost2 = sum(flip_cost(v,allele2) for v in variants[1])
		costs.append(cost1 + cost2)
	return min(costs)


def brute_force_phase(read_set, all_heterozygous):
	"""Solves MEC by enumerating all possible bipartitions."""
	assert len(read_set) < 10, "Too many reads for brute force"
	positions = read_set.get_positions()
	if all_heterozygous:
		possible_assignments = [(0,1), (1,0)]
	else:
		possible_assignments = [(0,0), (0,1), (1,0), (1,1)]
	# bit i in "partition" encodes to which set read i belongs
	best_partition = None
	best_cost = None
	solution_count = 0
	for partition in range(2**len(read_set)):
		print('Looking at partition {{:0>{}b}}'.format(len(read_set)).format(partition))
		# compute cost induced by that partition
		cost = 0
		for p in positions:
			# find variants covering this position
			variants = [[],[]]
			for n, read in enumerate(read_set):
				i = (partition >> n) & 1
				for variant in read:
					if variant.position == p:
						variants[i].append(variant)
			c = column_cost(variants, possible_assignments)
			print('    position: {}, variants: {} --> cost = {}'.format(p, str(variants), c))
			cost += c
		print('  --> cost for this partitioning:', cost)
		if (best_cost is None) or (cost < best_cost):
			best_partition = partition
			best_cost = cost
			solution_count = 1
		elif cost == best_cost:
			solution_count += 1
	# Each partition has its inverse with the same cost
	assert solution_count % 2 == 0
	return best_cost, [(best_partition>>x) & 1 for x in range(len(read_set))], solution_count//2


def compare_phasing(reads, all_heterozygous):
	"""Compares DPTable based phasing to brute force phasing and returns string representation of superreads."""
	rs = string_to_readset(reads)
	dp_table = DPTable(rs, all_heterozygous)
	superreads = dp_table.get_super_reads()
	assert len(superreads) == 2
	assert len(superreads[0]) == len(superreads[1])
	for v1, v2 in zip(*superreads):
		assert v1.position == v2.position
	result = tuple(sorted(''.join(str(v.allele) for v in sr) for sr in superreads))
	cost = dp_table.get_optimal_cost()
	partition = dp_table.get_optimal_partitioning()
	expected_cost, expected_partition, solution_count = brute_force_phase(rs, all_heterozygous)
	inverse_partition = [1-p for p in partition]
	print()
	print(superreads[0])
	print(superreads[1])
	print('Partition:', partition)
	print('Expected: ', expected_partition)
	print('Result:')
	print(result[0])
	print(result[1])
	print('Cost:', cost)
	assert (partition == expected_partition) or (inverse_partition == expected_partition)
	assert solution_count == 1
	assert cost == expected_cost
	# TODO: compute expected haplotypes based on expected_partition and compare to superreads
	return result


def test_phase1():
	reads = """
	 10
	 010
	 010
	"""
	s1, s2 = compare_phasing(reads, True)
	assert s1 == '010' and s2 == '101'


def test_phase2():
	reads = """
	  1  11010
	  00 00101
	  001 0101
	"""
	s1, s2 = compare_phasing(reads, True)
	assert s1 == '00100101' and s2 == '11011010'


def test_phase3():
	reads = """
	  1  11010
	  00 00101
	  001 01010
	"""
	s1, s2 = compare_phasing(reads, True)
	assert s1 == '001001010' and s2 == '110110101'


def test_phase4():
	reads = """
	  1  11010
	  00 00101
	  001 01110
	   1    111
	"""
	s1, s2 = compare_phasing(reads, True)

#TODO: test cases for weighted version (all base qualities are set to 1 right now)

#TODO: test cases to test "ties" (equals scores)

#TODO: test cases for all_heterozygous off
