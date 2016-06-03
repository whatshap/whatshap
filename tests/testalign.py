"""
Copied from https://bitbucket.org/marcelm/sqt/src/af255d54a21815cb9a3e0b279b431a320d4626bd/tests/testalign.py
"""
from whatshap.align import edit_distance as ed
from random import choice, seed, randint
from nose.tools import raises

STRING_PAIRS = [
	('', ''),
	('', 'A'),
	('A', 'A'),
	('AB', ''),
	('AB', 'ABC'),
	('TGAATCCC', 'CCTGAATC'),
	('ANANAS', 'BANANA'),
	('SISSI', 'MISSISSIPPI'),
	('GGAATCCC', 'TGAGGGATAAATATTTAGAATTTAGTAGTAGTGTT'),
	('TCTGTTCCCTCCCTGTCTCA', 'TTTTAGGAAATACGCC'),
	('TGAGACACGCAACATGGGAAAGGCAAGGCACACAGGGGATAGG', 'AATTTATTTTATTGTGATTTTTTGGAGGTTTGGAAGCCACTAAGCTATACTGAGACACGCAACAGGGGAAAGGCAAGGCACA'),
	('TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA', 'TTTTAGGAAATACGCCTGGTGGGGTTTGGAGTATAGTGAAAGATAGGTGAGTTGGTCGGGTG'),
	('A', 'TCTGCTCCTGGCCCATGATCGTATAACTTTCAAATTT'),
	('GCGCGGACT', 'TAAATCCTGG'),
	]


seed(10)

def randstring():
	return ''.join(choice('AC') for _ in range(randint(0, 10)))

STRING_PAIRS.extend((randstring(), randstring()) for _ in range(100000))


def test_edit_distance():
	assert ed('', '') == 0
	assert ed('', 'A') == 1
	assert ed('A', 'B') == 1
	assert ed('A', 'A') == 0
	assert ed('A', 'AB') == 1
	assert ed('BA', 'AB') == 2
	for s, t in STRING_PAIRS:
		assert ed(s, '') == len(s)
		assert ed('', s) == len(s)
		assert ed(s, t) == ed(t, s)


def test_edit_distance_bytes():
	assert ed(b'', b'') == 0
	assert ed(b'', b'A') == 1
	assert ed(b'A', b'B') == 1
	assert ed(b'A', b'A') == 0
	assert ed(b'A', b'AB') == 1
	assert ed(b'BA', b'AB') == 2
	for s, t in STRING_PAIRS:
		s = s.encode('ascii')
		t = t.encode('ascii')
		assert ed(s, b'') == len(s)
		assert ed(b'', s) == len(s)
		assert ed(s, t) == ed(t, s)


def assert_banded(s, t, maxdiff):
	banded_dist = ed(s, t, maxdiff=maxdiff)
	true_dist = ed(s, t)
	if true_dist > maxdiff:
		assert banded_dist > maxdiff
	else:
		assert banded_dist == true_dist


def test_edit_distance_banded():
	for maxdiff in range(5):
		assert_banded('ABC', '', maxdiff)

		for s, t in STRING_PAIRS:
			assert_banded(s, '', maxdiff)
			assert_banded('', s, maxdiff)
			assert_banded(s, t, maxdiff)
			assert_banded(t, s, maxdiff)
