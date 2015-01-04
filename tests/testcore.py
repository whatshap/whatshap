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
	r.addVariant(100, 'C', 0, 23)
	rs.add(r)

	r = Read('Read C', 17)
	r.addVariant(102, 'G', 1, 27)
	r.addVariant(105, 'T', 0, 14)
	rs.add(r)

	rs.finalize()
	assert len(rs) == 3


def test_phase_empty_readset():
	rs = ReadSet()
	rs.finalize()
	dp_table = DPTable(rs, all_heterozygous=False)
	superreads = dp_table.getSuperReads()
