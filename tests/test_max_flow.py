from phasingutils import string_to_readset
from whatshap.maxflow import Node,one_d_range_tree

def test_selection():
	reads = string_to_readset("""
	  11
	  000
	    11
	      00
	    11
	""")
	one_d_range_tree(reads)
	assert 0==1


#def test_selection2():
#	reads = string_to_readset("""
#	  1111
#	     111
#	     1  111
#	     1     11
#	    1      11
#	""")
