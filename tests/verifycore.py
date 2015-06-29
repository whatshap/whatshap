import sys
from nose.tools import raises
from whatshap.core import Read, DPTable, ReadSet, Variant
from whatshap.verification import verify_mec_score_and_partitioning
from .phasingutils import string_to_readset, matrix_to_readset, brute_force_phase


def verify(rs, all_heterozygous = False) :
	dp_table = DPTable(rs, all_heterozygous)
	verify_mec_score_and_partitioning(dp_table, rs)


def test_string() :
	reads = """
	  0             0
	  110111111111
	  00100
	       0001000000
	       000
	        10100
	              101
	"""
	rs = string_to_readset(reads)
	verify(rs, True)
	verify(rs, False)


def test_matrix() :
	rs = matrix_to_readset(open('tests/test.matrix','r'))
	verify(rs, True)
	verify(rs, False)
