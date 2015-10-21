from phasingutils import string_to_readset
from whatshap.maxflow import Leaf_node,one_d_range_tree

def test_selection():
	reads = string_to_readset("""
	  11
	  000
	    11
	      00
	    11
	""")
	first_tree_attemp=one_d_range_tree(reads)
	#first_tree_attemp.build_list(reads)
	assert 0==1


#def test_selection2():
#	reads = string_to_readset("""
#	  1111
#	     111
#	     1  111
#	     1     11
#	    1      11
#	""")
