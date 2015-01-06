from whatshap.core import Read, DPTable, ReadSet

def test_read():
	r = Read("name", 15)
	assert r.getName() == "name"
	assert r.getMapqs()[0] == 15


def test_empty_readset():
	rs = ReadSet()
	assert len(rs) == 0


def test_readset():
	rs = ReadSet()
	r = Read('Read A', 56)
	r.addVariant(100, 'A', 1, 37)
	r.addVariant(101, 'C', 0, 18)
	rs.add(r)

	r = Read('Read B', 0)
	r.addVariant(101, 'C', 0, 23)
	rs.add(r)

	r = Read('Read C', 17)
	r.addVariant(99, 'G', 1, 27)
	r.addVariant(105, 'T', 0, 14)
	rs.add(r)

	assert rs[0].getName() == 'Read A'
	assert rs[1].getName() == 'Read B'
	assert rs[2].getName() == 'Read C'

	rs.finalize()

	# should be sorted after finalization
	assert rs[0].getName() == 'Read C'
	assert rs[1].getName() == 'Read A'
	assert rs[2].getName() == 'Read B'

	assert len(rs) == 3

	assert rs.getPositions() == [99,100,101,105]

	r = rs.getByName('Read A')
	assert r.getName() == 'Read A'
	assert r.getMapqs() == (56,), str(r.getMapqs())

	r = rs.getByName('Read B')
	assert r.getName() == 'Read B'
	assert r.getMapqs() == (0,)

	try:
		# Should raise a KeyError for non-existing read name
		r = rs.getByName('foo')
		assert False
	except KeyError:
		pass

	# TODO: Test subset method



def test_phase_empty_readset():
	rs = ReadSet()
	rs.finalize()
	dp_table = DPTable(rs, all_heterozygous=False)
	superreads = dp_table.getSuperReads()
