from nose.tools import raises
from whatshap.core import Read, DPTable, ReadSet, Variant
from phasingutils import string_to_readset, brute_force_phase
from whatshap.verification import verify_mec_score_and_partitioning


def verify(reads, all_heterozygous = False) :
	rs = string_to_readset(reads, None)
	dp_table = DPTable(rs, all_heterozygous)
	super_reads = dp_table.get_super_reads()
	verify_mec_score_and_partitioning(dp_table, super_reads)


def test_phase_trivial() :
	reads = """
          1 1
          000
           1
        """
	verify(reads, True)
	verify(reads, False)
