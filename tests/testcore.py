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


def phase_from_string(s):
	"""
	Phase the reads specified via the string s and return the result as a
	pair of strings.
	"""
	s = textwrap.dedent(s).strip()
	bits = { '0': 'A', '1': 'C', 'E': 'G' }
	rs = ReadSet()
	for index, line in enumerate(s.split('\n'), 1):
		read = Read('Read {}'.format(index), 50)
		for pos, c in enumerate(line, 1):
			if c == ' ':
				continue
			read.add_variant(position=pos * 10, base=bits[c], allele=int(c), quality=22)
		rs.add(read)

	print(rs)
	dp_table = DPTable(rs, all_heterozygous=True)
	superreads = dp_table.get_super_reads()
	assert len(superreads) == 2
	assert len(superreads[0]) == len(superreads[1])
	for v1, v2 in zip(*superreads):
		assert v1.position == v2.position
	print()
	print(superreads[0])
	print(superreads[1])
	result = tuple(''.join(str(v.allele) for v in sr) for sr in superreads)
	print('Result:')
	print(result[0])
	print(result[1])
	return result


def test_phase():
	reads = """
	 10
	 010
	 010
	"""
	s1, s2 = phase_from_string(reads)
	assert (s1 == '101' and s2 == '010') or (s1 == '010' and s2 == '101')


# TODO
def test_phase2():
	IN = """
	 0011  1101
	  1100 00101
	   1001 01010
	"""
	IN = """
	 0011  11010
	  1100 00101
	   1001001010
	"""
