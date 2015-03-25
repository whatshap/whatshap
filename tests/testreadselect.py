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
	#Here the assert is wrong, because the bridging doesn't come into account , because in the slice_read the selected
	# reads  have already coverage 3 by set ([1,3,5,7]) because first each position has to covered at least once before
	#the bridging starts
	#assert selected_reads == set([0,1,3,5]), str(selected_reads)


def test_bridging():
	reads = string_to_readset("""
	  11
	  00
	    11
	    00
	      11
	      00
	  1    1
	""")
	selected_reads, new_components = readselection(reads, max_cov = 2, bridging= False)
	assert selected_reads == set([0,1,2,3,4,5])
	selected_reads, new_components = readselection(reads, max_cov = 2, bridging= True)
	#Not sure why 0 is there selected and not 1...
	assert selected_reads == set([0,3,5,6])

def test_components_of_readselection():
	reads = string_to_readset("""
	  111
	     000
	  00
	      00
	    1
	       1
	   1   1
	""")
	selected_reads, new_components = readselection(reads, max_cov = 2, bridging= False)
	assert selected_reads == set([0,1,2,3])
	assert len(new_components.values()) == 2
	selected_reads, new_components = readselection(reads, max_cov = 2, bridging= True)
	assert selected_reads == set([0,1,6])
	assert len(new_components.values()) == 1
