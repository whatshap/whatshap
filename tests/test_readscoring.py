"""
Test ReadScoring
"""

from whatshap.readscoring import score
from whatshap.core import Read, ReadSet, Variant

def test_readscoring_toy():
	readset = ReadSet()
	read1 = Read('name1', 15)
	read1.add_biallelic_variant(0, 0, 1)
	read1.add_biallelic_variant(1, 0, 1)
	read1.add_biallelic_variant(2, 0, 1)
	read1.add_biallelic_variant(3, 1, 1)
	readset.add(read1)
	read2 = Read('name2', 15)
	read2.add_biallelic_variant(1, 1, 1)
	read2.add_biallelic_variant(2, 0, 1)
	read2.add_biallelic_variant(3, 0, 1)
	read2.add_biallelic_variant(4, 1, 1)
	readset.add(read2)
	read3 = Read('name3', 15)
	read3.add_biallelic_variant(2, 0, 1)
	read3.add_biallelic_variant(3, 1, 1)
	read3.add_biallelic_variant(4, 0, 1)
	read3.add_biallelic_variant(5, 1, 1)
	readset.add(read3)
	read4 = Read('name4', 15)
	read4.add_biallelic_variant(3, 0, 1)
	read4.add_biallelic_variant(4, 1, 1)
	read4.add_biallelic_variant(5, 0, 1)
	read4.add_biallelic_variant(6, 0, 1)
	readset.add(read4)
	read5 = Read('name5', 15)
	read5.add_biallelic_variant(4, 0, 1)
	read5.add_biallelic_variant(5, 1, 1)
	read5.add_biallelic_variant(6, 1, 1)
	read5.add_biallelic_variant(7, 0, 1)
	readset.add(read5)
	read6 = Read('name6', 15)
	read6.add_biallelic_variant(5, 0, 1)
	read6.add_biallelic_variant(6, 0, 1)
	read6.add_biallelic_variant(7, 0, 1)
	read6.add_biallelic_variant(8, 1, 1)
	readset.add(read6)
	read7 = Read('name7', 15)
	read7.add_biallelic_variant(6, 1, 1)
	read7.add_biallelic_variant(7, 0, 1)
	read7.add_biallelic_variant(8, 0, 1)
	read7.add_biallelic_variant(9, 1, 1)
	readset.add(read7)
	sim = score(readset, 2, 0.05, 1)

	assert sim[0][0] < 0.0
	assert sim[0][1] > 0.0
	assert sim[0][2] < 0.0
	assert sim[0][3] > 0.0
	assert sim[0][4] < 0.0
	assert sim[0][5] > 0.0
	assert sim[1][0] < 0.0
	assert sim[1][1] > 0.0
	assert sim[1][2] < 0.0
	assert sim[1][3] > 0.0
	assert sim[1][4] < 0.0
	assert sim[2][0] < 0.0
	assert sim[2][1] > 0.0
	assert sim[2][2] < 0.0
	assert sim[2][3] > 0.0
	assert sim[3][0] < 0.0
	assert sim[3][1] > 0.0
	assert sim[3][2] < 0.0
	assert sim[4][0] < 0.0
	assert sim[4][1] > 0.0
	assert sim[5][0] < 0.0
