"""
Test Read and ReadSet classes
"""
from nose.tools import raises
from whatshap.core import Read, ReadSet, Variant


def test_read():
	r = Read("name", 15)
	assert r.name == "name"
	assert r.mapqs[0] == 15

	assert r.is_sorted()

	r.add_variant(100, 1, 37)
	r.add_variant(23, 0, 99)
	assert not r.is_sorted()
	r.sort()
	assert r.is_sorted()

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
