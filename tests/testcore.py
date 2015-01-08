from whatshap.core import Read, DPTable, ReadSet

def test_read():
	r = Read("name", 15)
	assert r.get_name() == "name"
	assert r.get_mapqs()[0] == 15

	assert r.is_sorted()

	r.add_variant(100, 'A', 1, 37)
	r.add_variant(23, 'T', 0, 99)
	assert not r.is_sorted()
	r.sort()
	assert r.is_sorted()

	assert 100 in r
	assert 23 in r
	assert not 22 in r
	assert not 24 in r
	assert not 1000 in r
	assert not -1000 in r


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

	assert rs[0].get_name() == 'Read A'
	assert rs[1].get_name() == 'Read B'
	assert rs[2].get_name() == 'Read C'

	rs.sort()

	# should be sorted after finalization
	assert rs[0].get_name() == 'Read C'
	assert rs[1].get_name() == 'Read A'
	assert rs[2].get_name() == 'Read B'

	assert len(rs) == 3

	assert rs.get_positions() == [99, 100, 101, 105]

	r = rs['Read A']
	assert r.get_name() == 'Read A'
	assert r.get_mapqs() == (56,), str(r.get_mapqs())

	r = rs['Read B']
	assert r.get_name() == 'Read B'
	assert r.get_mapqs() == (0,)

	try:
		# Should raise a KeyError for non-existing read name
		r = rs['foo']
		assert False
	except KeyError:
		pass

	# TODO: Test subset method



def test_phase_empty_readset():
	rs = ReadSet()
	dp_table = DPTable(rs, all_heterozygous=False)
	superreads = dp_table.get_super_reads()
