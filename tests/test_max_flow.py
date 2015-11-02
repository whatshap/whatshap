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


def test_algo_by_mini_readset():
	reads = string_to_readset("""
	  11
	  00
	    11
	""")
	first_tree_attemp=one_d_range_tree(reads)
	#complete tree is the root node and from there on it gets down to the leafs
	complete_tree=first_tree_attemp.get_complete_tree()
	root_coverage=complete_tree.get_coverage()
	assert root_coverage==(1,2)

	#assert 0==1

#def test_struct_of_tree():
#	reads = string_to_readset("""
#	  11
#	  00
#	    11
#	""")
#	first_tree_attemp=one_d_range_tree(reads)
#	root_node=first_tree_attemp.get_complete_tree()
#	new_root_node=root_node.return_node()
#	print('Find out if we get the leaf nodes of the root node')
	#Method_list=new_root_node.get_all_leaf_nodes_of_subtree()



#	first_cov=root_node.get_coverage()
#	print(first_cov)
#	leaf_list=first_tree_attemp.get_leaf_list_of_tree()
#	print('Leaf List ')
#	print(len(leaf_list))
#	print(leaf_list)
#	assert  (len(leaf_list)==4)
#	for i,r in enumerate(reads):
#		print('I of enumerate %d '%i)
#		print(leaf_list[i].get_parent().get_coverage())
#		one_leaf=leaf_list[i].get_sibling()
#		print('one_leaf')
#		print(one_leaf)
#		firs_pos=r[0]
#		print('First position of the read')
#		print(firs_pos.position)
		#Leaf_node_of_first_position=leaf_list.index(int(firs_pos.position))
		#last_pos=read[len(read)-1]


#	for i,leaf in enumerate(leaf_list):
#		l_parent=leaf.get_parent()
#		r_child=l_parent.get_right_child()

#	assert 0==1



def test_split_node():
	reads = string_to_readset("""
	  11
	  000
	   11
	    000
	     111
	""")
	tree=one_d_range_tree(reads)
	root_node=tree.get_complete_tree()
	leaf_list_of_tree=tree.get_leaf_list_of_tree()
	for i,l_node in enumerate(leaf_list_of_tree):
		print('Index of leaf_list %d'%i)
		print(l_node.get_coverage())
		print(l_node.get_parent().get_coverage())

	assert 0==1
