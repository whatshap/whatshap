"""
Edit distance tests

The initial version was copied from https://bitbucket.org/marcelm/sqt/src/af255d54a21815cb9a3e0b279b431a320d4626bd/tests/testalign.py
"""
from whatshap.align import edit_distance as ed
from whatshap.align import edit_distance_affine_gap as ed_aff
from random import choice, seed, randint

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


def randstring():
	return ''.join(choice('AC') for _ in range(randint(0, 10)))


seed(10)
STRING_PAIRS.extend((randstring(), randstring()) for _ in range(1000))


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


def test_edit_distance_affine():
	for mismatch_cost in [1, 10, 30, 40, 50]:
		for gap_start in [1, 10, 30, 40, 50]:
			assert ed_aff('', '', [], gap_start, 10) == 0
			assert ed_aff('', 'A', [], gap_start, 10) == gap_start
			assert ed_aff('A', 'B', [mismatch_cost], gap_start, 10) == min(gap_start * 2,
				mismatch_cost)
			assert ed_aff('A', 'A', [mismatch_cost], gap_start, 10) == 0
			assert ed_aff('A', 'AB', [mismatch_cost], gap_start, 10) == gap_start
			assert ed_aff('BA', 'AB', [mismatch_cost] * 2, gap_start, 100) == min(2 * mismatch_cost,
				2 * gap_start)
			for s, t in STRING_PAIRS:
				if s != '':
					assert ed_aff(s, '', [mismatch_cost] * len(s), gap_start, 10) == gap_start + (
								len(s) - 1) * 10
					assert ed_aff('', s, [], gap_start, 10) == gap_start + (len(s) - 1) * 10
				assert ed_aff(s, t, [mismatch_cost] * len(s), gap_start, 10) == ed_aff(t, s,
					[mismatch_cost] * len(t), gap_start, 10)


def test_edit_distance_affine_bytes():
	for mismatch_cost in [1, 10, 20, 30, 40, 50]:
		for gap_start in [1, 10, 20, 30, 40, 50]:
			assert ed_aff(b'', b'', [], gap_start, 10) == 0
			assert ed_aff(b'', b'A', [], gap_start, 10) == gap_start
			assert ed_aff(b'A', b'B', [mismatch_cost], gap_start, 10) == min(gap_start * 2,
				mismatch_cost)
			assert ed_aff(b'A', b'A', [mismatch_cost], gap_start, 10) == 0
			assert ed_aff(b'A', b'AB', [mismatch_cost], gap_start, 10) == gap_start
			assert ed_aff(b'BA', b'AB', [mismatch_cost] * 2, gap_start, 100) == min(
				2 * mismatch_cost, 2 * gap_start)
			for s, t in STRING_PAIRS:
				s = s.encode('ascii')
				t = t.encode('ascii')
				if s != b'':
					assert ed_aff(s, b'', [mismatch_cost] * len(s), gap_start, 10) == gap_start + (
								len(s) - 1) * 10
					assert ed_aff(b'', s, [], gap_start, 10) == gap_start + (len(s) - 1) * 10
				assert ed_aff(s, t, [mismatch_cost] * len(s), gap_start, 10) == ed_aff(t, s,
					[mismatch_cost] * len(t), gap_start, 10)


def test_mismatches():
	for i in range(10):
		rand_costs = [randint(10, 70) for j in range(5)]
		assert ed_aff('AAAAA', 'TTTTT', rand_costs, 100, 100) == sum(rand_costs)
		assert ed_aff('ATGCT', 'ATCCT', rand_costs, 100, 100) == rand_costs[2]
		assert ed_aff('ATGGA', 'ATGTTCA', rand_costs, 80, 10) == rand_costs[3] + 80 + 10


def test_small_examples():
	assert ed_aff('AGTCCGGTG', 'AGTCCATCGGTC', [30, 40, 20, 20, 50, 60, 10, 20, 5], 40, 10) == 65
	assert ed_aff('ATGGCCG', 'ATCGCTG', [40, 50, 10, 40, 50, 10, 40], 20, 10) == 20
	assert ed_aff('ATCCTC', 'ATCGGGCTC', [50] * 6, 10, 5) == 20


def test_compare_to_edit_dist():
	for s, t in STRING_PAIRS:
		assert ed(s, t) == ed_aff(s, t, [1] * len(s), 1, 1)
