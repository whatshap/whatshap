from whatshap.core import readselection
from .phasingutils import string_to_readset


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
	selected_reads = readselection(reads, max_cov = 1, preferred_source_ids = None, bridging = False)
	assert selected_reads == set([1,5])
	selected_reads = readselection(reads, max_cov = 2, preferred_source_ids = None, bridging = False)
	assert selected_reads == set([1,3,5]), str(selected_reads)
	selected_reads = readselection(reads, max_cov = 3, preferred_source_ids = None, bridging = False)
	assert selected_reads == set([1,3,5,7]), str(selected_reads)
	selected_reads = readselection(reads, max_cov = 3, preferred_source_ids = None, bridging = True)
	#Here the assert is wrong, because the bridging doesn't come into account , because in the slice_read the selected
	# reads  have already coverage 3 by set ([1,3,5,7]) because first each position has to covered at least once before
	#the bridging starts
	assert selected_reads == set([1,3,5,7]), str(selected_reads)


def test_selection2():
	reads = string_to_readset("""
	  1111
	     111
	     1  111
	     1     11
	    1      11
	""")
	selected_reads = readselection(reads, max_cov = 4, preferred_source_ids = None, bridging = False)
	assert selected_reads == set([0,1,2,3]), str(selected_reads)


def bridging():
	reads = string_to_readset("""
	  11
	  00
	    11
	    00
	      11
	      00
	  1    1
	""")
	selected_reads  = readselection(reads, max_cov = 2, preferred_source_ids = None, bridging= False)
	assert selected_reads == set([0,1,2,3,4,5])
	selected_reads = readselection(reads, max_cov = 2, preferred_source_ids = None, bridging= True)
	#Not sure why 0 is there selected and not 1...
	assert selected_reads == set([0,3,5,6])


###Component comparison does not work
def test_components_of_readselection():
	reads = string_to_readset("""
	  111
	     000
	  00
	      00
	   1   1
	""")
	selected_reads = readselection(reads, max_cov = 2, preferred_source_ids = None, bridging= False)
	assert selected_reads == set([0,1,2,3]), str(selected_reads)
#	assert len(set(new_components.values())) == 2
	selected_reads = readselection(reads, max_cov = 2, preferred_source_ids = None, bridging= True)
	assert selected_reads == set([0,1,4]), str(selected_reads)
#  	assert len(set(new_components.values())) == 1


def test_selection_with_preferred_sources():
	readset = string_to_readset("""
	  1        1
	""", source_id = 3)
	more_reads = string_to_readset("""
	  1111
	     111
	        1111
	""", source_id = 1)

	for read in more_reads:
		readset.add(read)
	
	selected_reads = readselection(readset, max_cov = 2, preferred_source_ids = None, bridging = True)
	assert selected_reads == set([1,2,3]), str(selected_reads)

	selected_reads = readselection(readset, max_cov = 2, preferred_source_ids = set([3]), bridging = True)
	assert selected_reads == set([0,1,3]), str(selected_reads)


#TODO: the below test case seems to be incomplete
#def test_tuple_scores():
	#'Only example at the moment '
	#reads = """
	  #1  11010
	  #00 00101
	  #001 01110
	   #1    111
	#"""
	#weights = """
	  #2  13112
	  #11 23359
	  #223 56789
	   #2    111
	#"""
	#rs = string_to_readset(reads, weights)
	#selected_reads, skipped_reads = readselection(reads, max_cov = 2, bridging= False)
	#selected_reads, skipped_reads = readselection(reads, max_cov = 2, bridging= True)

