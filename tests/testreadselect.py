from phasingutils import string_to_readset
from whatshap.readselect import readselection

def test_selection():
	reads = string_to_readset("""
	  1  1
	  00
	  0   1
	  10  1
	  1   1
	    11
	  0   1
	  1    1
	""")
	selected_reads, new_components = readselection(reads, max_cov = 1, bridging = False)
	assert selected_reads == set([1,5])
	selected_reads, new_components = readselection(reads, max_cov = 2, bridging = False)
	assert selected_reads == set([1,3,5]), str(selected_reads)
	selected_reads, new_components = readselection(reads, max_cov = 3, bridging = False)
	assert selected_reads == set([1,3,5,7]), str(selected_reads)
	selected_reads, new_components = readselection(reads, max_cov = 3, bridging = True)
	assert selected_reads == set([0,1,3,5]), str(selected_reads)

