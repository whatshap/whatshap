from whatshap.core import Read, DPTable, ReadSet

def test_read():
	r = Read("name", 15)
	assert r.getName() == "name"
	assert r.getMapqs()[0] == 15


def test_readset():
	rs = ReadSet()
	assert len(rs) == 0


def test_empty_readset():
	rs = ReadSet()
	rs.finalize()
	dp_table = DPTable(rs, all_heterozygous=False)
	# TODO C++ exception in the following line
	superreads = dp_table.getSuperReads()
