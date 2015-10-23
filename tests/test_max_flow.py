from phasingutils import string_to_readset
from whatshap.maxflow import Leaf_node,one_d_range_tree

def test_tree_struct():
	reads = string_to_readset("""
	  11
	  000
	    11
	      00
	    11
	""")
	first_tree_attemp=one_d_range_tree(reads)
	#complete tree is the root node and from there on it gets down to the leafs
	complete_tree=first_tree_attemp.get_complete_tree()
	root_coverage=complete_tree.get_coverage()
	assert root_coverage==(1,3)
	left_child=complete_tree.get_left_child()
	right_child=complete_tree.get_right_child()
	assert left_child.get_coverage()==(2,3)
	assert right_child.get_coverage()==(1,2)
	assert left_child.get_balance()==0
	assert right_child.get_balance()==0
	assert left_child.get_min_coverage()==2
	assert right_child.get_max_coverage()==2

	#assert 0==1


#def test_selection2():
#	reads = string_to_readset("""
#	  1111
#	     111
#	     1  111
#	     1     11
#	    1      11
#	""")
