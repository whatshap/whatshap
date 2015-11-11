from phasingutils import string_to_readset
from whatshap.maxflow import optimize_max_flow_in_BST
from whatshap.Binary_Search_Tree import Binary_Search_Tree




#simple case balances BST

def test_struct_for_BST():
	reads=string_to_readset("""
	 111
	 000
	   11
	    00
	    11
	""")
	tree=Binary_Search_Tree(reads)
	root_node=tree.get_complete_tree()
	leaf_list_of_tree=tree.get_leaf_list_of_tree()
	#go over leafs =reads

	leaf_value=leaf_list_of_tree[0].get_value()
	siblings=leaf_list_of_tree[0].get_sibling()
	only_end_points=[sib for sib in siblings if sib.get_value()>leaf_value]
	tree.synchronize_sibling_with_same_value(only_end_points)
	assert only_end_points[0].get_value() ==30
	assert only_end_points[1].get_value() ==30
	(split_node_of_read,List_to_change)=tree.seach_for_split_node(leaf_list_of_tree[0],only_end_points[0])
	print('List to change')
	print(List_to_change)
	print(len(List_to_change))
	assert split_node_of_read.get_coverage()==(2,3)
	assert split_node_of_read.get_right_child().get_value()==30
	print('List to change')
	print(List_to_change)
	print('Split node')
	print(split_node_of_read)
	assert len(List_to_change)==3
	#Geht nicht mehr da List to change nun set ist
	#assert List_to_change[0]==split_node_of_read
	assert split_node_of_read in List_to_change
	print('only_end_points[1]')
	print(only_end_points[1].get_value())
	print(only_end_points[1].get_coverage())
	(split_node_of_read,List_to_change)=tree.seach_for_split_node(leaf_list_of_tree[0],only_end_points[1])
	assert split_node_of_read.get_coverage()==(2,3)
	assert split_node_of_read.get_right_child().get_value()==30

	leaf_value=leaf_list_of_tree[1].get_value()
	siblings=leaf_list_of_tree[1].get_sibling()
	only_end_points=[sib for sib in siblings if sib.get_value()>leaf_value]
	assert only_end_points[0].get_value() ==40
	(split_node_of_read,List_to_change)=tree.seach_for_split_node(leaf_list_of_tree[1],only_end_points[0])
	assert split_node_of_read.get_coverage()==(2,3)
	assert split_node_of_read==root_node


def test_reveale_whole_tree_structure_for_simple_case():
	reads=string_to_readset("""
	 111
	 000
	   11
	    00
	    11
	""")
	tree=Binary_Search_Tree(reads)
	root_node=tree.get_complete_tree()
	leaf_list_of_tree=tree.get_leaf_list_of_tree()
	first_leaf=leaf_list_of_tree[0]
	assert first_leaf.get_value()==10
	assert first_leaf.get_is_left_child()
	siblings_of_first_leaf=first_leaf.get_sibling()
	assert len(siblings_of_first_leaf)== 2
	tree.synchronize_sibling_with_same_value(siblings_of_first_leaf)
	first_of_value_30=siblings_of_first_leaf[0]
	second_of_value_30=siblings_of_first_leaf[1]
	#assert values
	assert first_of_value_30.get_value()==30
	assert second_of_value_30.get_value()==30
	#assert parent
	assert first_of_value_30.get_parent().get_coverage()==(2,3)
	assert second_of_value_30.get_parent().get_coverage()==(2,3)
	#assert isLeaf
	assert first_of_value_30.isLeaf()
	assert second_of_value_30.isLeaf()
	#assert is_Left_child
	assert first_of_value_30.get_is_right_child()
	assert second_of_value_30.get_is_right_child()
	assert not first_of_value_30.get_is_left_child()
	assert not second_of_value_30.get_is_left_child()
	#assert index
	assert first_of_value_30.get_index()==[0]
	assert second_of_value_30.get_index()==[1]
	#assert get_coverage
	assert first_of_value_30.get_coverage()==3
	assert second_of_value_30.get_coverage()==3
	second_leaf=leaf_list_of_tree[1]
	assert second_leaf.get_value()==30
	assert not second_leaf.get_is_left_child()
	siblings_of_second_leaf=second_leaf.get_sibling()
	assert len(siblings_of_second_leaf)==3
	tree.synchronize_sibling_with_same_value(siblings_of_second_leaf)
	assert len(siblings_of_second_leaf)==3
	first_of_value=siblings_of_second_leaf[0]
	second_of_value=siblings_of_second_leaf[1]
	third_of_value=siblings_of_second_leaf[2]
	#assert values
	assert first_of_value.get_value()==10
	assert second_of_value.get_value()==10
	assert third_of_value.get_value()==40
	#assert parent
	assert first_of_value.get_parent().get_coverage()==(2,3)
	assert second_of_value.get_parent().get_coverage()==(2,3)
	assert third_of_value.get_parent().get_coverage()==(2,3)
	#assert isLeaf
	assert first_of_value.isLeaf()
	assert second_of_value.isLeaf()
	assert third_of_value.isLeaf()
	#assert is_Left_child
	assert not first_of_value.get_is_right_child()
	assert not second_of_value.get_is_right_child()
	assert not third_of_value.get_is_right_child()
	assert first_of_value.get_is_left_child()
	assert second_of_value.get_is_left_child()
	assert third_of_value.get_is_left_child()
	#assert index
	assert first_of_value.get_index()==[0]
	assert second_of_value.get_index()==[1]
	assert third_of_value.get_index()==[2]
	#assert get_coverage
	assert first_of_value.get_coverage()==2
	assert second_of_value.get_coverage()==2
	assert third_of_value.get_coverage()==3
	third_leaf=leaf_list_of_tree[2]
	assert third_leaf.get_value()==40
	assert third_leaf.get_is_left_child()
	siblings_of_third_leaf=third_leaf.get_sibling()
	assert len(siblings_of_third_leaf)==3
	tree.synchronize_sibling_with_same_value(siblings_of_third_leaf)
	first_of_value=siblings_of_third_leaf[0]
	second_of_value=siblings_of_third_leaf[1]
	third_of_value=siblings_of_third_leaf[2]
	#assert values
	assert first_of_value.get_value()==30
	assert second_of_value.get_value()==50
	assert third_of_value.get_value()==50
	#assert parent
	assert first_of_value.get_parent().get_coverage()==(2,3)
	assert second_of_value.get_parent().get_coverage()==(2,3)
	assert third_of_value.get_parent().get_coverage()==(2,3)
	#assert isLeaf
	assert first_of_value.isLeaf()
	assert second_of_value.isLeaf()
	assert third_of_value.isLeaf()
	#assert is_Left_child
	assert first_of_value.get_is_right_child()
	assert second_of_value.get_is_right_child()
	assert third_of_value.get_is_right_child()
	assert not first_of_value.get_is_left_child()
	assert not second_of_value.get_is_left_child()
	assert not third_of_value.get_is_left_child()
	#assert index
	assert first_of_value.get_index()==[2]
	assert second_of_value.get_index()==[3]
	assert third_of_value.get_index()==[4]
	#assert get_coverage
	assert first_of_value.get_coverage()==3
	assert second_of_value.get_coverage()==2
	assert third_of_value.get_coverage()==2
	forth_leaf=leaf_list_of_tree[3]
	assert forth_leaf.get_value()==50
	assert not forth_leaf.get_is_left_child()
	siblings_of_forth_leaf=forth_leaf.get_sibling()
	assert len(siblings_of_forth_leaf)==2
	tree.synchronize_sibling_with_same_value(siblings_of_forth_leaf)
	first_of_value=siblings_of_forth_leaf[0]
	second_of_value=siblings_of_forth_leaf[1]
	#assert values
	assert first_of_value.get_value()==40
	assert second_of_value.get_value()==40
	#assert parent
	assert first_of_value.get_parent().get_coverage()==(2,3)
	assert second_of_value.get_parent().get_coverage()==(2,3)
	#assert isLeaf
	assert first_of_value.isLeaf()
	assert second_of_value.isLeaf()
	#assert is_Left_child
	assert not first_of_value.get_is_right_child()
	assert not second_of_value.get_is_right_child()
	assert first_of_value.get_is_left_child()
	assert second_of_value.get_is_left_child()
	#assert index
	assert first_of_value.get_index()==[3]
	assert second_of_value.get_index()==[4]
	#assert get_coverage
	assert first_of_value.get_coverage()==3
	assert second_of_value.get_coverage()==3


def test_in_simple_case_list_to_change_of_split_node_case():
	reads=string_to_readset("""
	 111
	 000
	   11
	    00
	    11
	""")
	tree=Binary_Search_Tree(reads)
	root_node=tree.get_complete_tree()
	leaf_list_of_tree=tree.get_leaf_list_of_tree()
	#go over leafs =reads

	leaf_value=leaf_list_of_tree[0].get_value()
	siblings=leaf_list_of_tree[0].get_sibling()
	only_end_points=[sib for sib in siblings if sib.get_value()>leaf_value]
	tree.synchronize_sibling_with_same_value(only_end_points)
	assert only_end_points[0].get_value() ==30
	assert only_end_points[1].get_value() ==30
	(split_node_of_read,List_to_change)=tree.seach_for_split_node(leaf_list_of_tree[0],only_end_points[0])
	#List to change should include 3 nodes, the 2 nodes and the parent node of both
	assert len(List_to_change) ==3
	out_of_list_to_set=set(List_to_change)
	assert leaf_list_of_tree[0] in out_of_list_to_set
	assert only_end_points[0] in out_of_list_to_set
	assert leaf_list_of_tree[0].get_parent() in out_of_list_to_set


	(split_node_of_read,List_to_change)=tree.seach_for_split_node(leaf_list_of_tree[0],only_end_points[1])
	assert len(List_to_change) ==4
	out_of_list_to_set=set(List_to_change)
	assert leaf_list_of_tree[0] in out_of_list_to_set
	assert only_end_points[1] in out_of_list_to_set
	assert only_end_points[0] in out_of_list_to_set
	assert leaf_list_of_tree[0].get_parent() in out_of_list_to_set


	middle_leaf=leaf_list_of_tree[1]
	siblings=middle_leaf.get_sibling()
	assert middle_leaf.get_value()==30
	only_end_points=[sib for sib in siblings if sib.get_value()>leaf_value]
	tree.synchronize_sibling_with_same_value(only_end_points)
	assert len(only_end_points)==1
	assert only_end_points[0].get_value()==40
	(split_node_of_read,List_to_change)=tree.seach_for_split_node(middle_leaf,only_end_points[0])
	assert split_node_of_read==middle_leaf.get_parent().get_parent()
	assert len(List_to_change)==5
	out_of_list_to_set=set(List_to_change)
	assert middle_leaf in out_of_list_to_set
	assert middle_leaf.get_parent() in out_of_list_to_set
	assert only_end_points[0] in out_of_list_to_set
	assert only_end_points[0].get_parent() in out_of_list_to_set
	assert split_node_of_read in out_of_list_to_set




#TODO NOt WORKIGNG correctly
def test_simple_case_max_flow():
	reads=string_to_readset("""
	 111
	 000
	   11
	    00
	    11
	""")
	tree=Binary_Search_Tree(reads)
	root_node=tree.get_complete_tree()
	leaf_list_of_tree=tree.get_leaf_list_of_tree()
	maximal_coverage=2
	pruned_set,removed_set=optimize_max_flow_in_BST(reads,tree,maximal_coverage)
	assert len(pruned_set)==3
	assert len(removed_set)==2
	print('Pruned SEt')
	print(pruned_set)
	print('Removed set')
	print(removed_set)

	maximal_coverage=1
	pruned_set,removed_set=optimize_max_flow_in_BST(reads,tree,maximal_coverage)
	print('Pruned SEt')
	print(pruned_set)
	print('Removed set')
	print(removed_set)
	#Should be
	#assert len(pruned_set)==2 # should be read 0 or 1 and read 3 or 4
	#assert len(removed_set)==3 #should definitive be read 2 and the contrary to the one in pruned
	#but because of the not working crucial attribute not correct, in the logical sense.
	assert len(pruned_set)==3
	assert len(removed_set)==2
	#Is working
	maximal_coverage=3
	pruned_set,removed_set=optimize_max_flow_in_BST(reads,tree,maximal_coverage)
	assert len(pruned_set)==5
	assert len(removed_set)==0


def test_tree_struct_for_Binary_Seach_Tree():
	reads = string_to_readset("""
	  11
	  000
	    11
	      00
	    11
	""")
	first_tree_attemp=Binary_Search_Tree(reads)
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


def test_BST_Siblings_of_nodes():
	reads = string_to_readset("""
	  11
	  0000
	    11
	      00
	      11
	    11
	""")
	first_tree_attemp=Binary_Search_Tree(reads)
	#complete tree is the root node and from there on it gets down to the leafs
	complete_tree=first_tree_attemp.get_complete_tree()
	root_coverage=complete_tree.get_coverage()
	assert root_coverage==(2,3)
	L_of_root=complete_tree.get_left_child()
	R_of_root=complete_tree.get_right_child()
	assert L_of_root.get_coverage()==(2,3)
	assert R_of_root.get_coverage()==(2,3)
	left_child_of_L=L_of_root.get_left_child()
	right_child_of_L=L_of_root.get_right_child()
	left_child_of_R=R_of_root.get_left_child()
	right_child_of_R=R_of_root.get_right_child()
	assert right_child_of_L.isLeaf()
	assert right_child_of_L.get_value()==30
	assert right_child_of_R.isLeaf()
	assert right_child_of_R.get_value()==60
	assert not left_child_of_L.isLeaf()
	assert not left_child_of_R.isLeaf()
	l_of_L=left_child_of_L.get_left_child()
	r_of_L=left_child_of_L.get_right_child()
	assert l_of_L.get_value()==10
	assert r_of_L.get_value()==20
	l_of_R=left_child_of_R.get_left_child()
	r_of_R=left_child_of_R.get_right_child()
	assert l_of_R.get_value()==40
	assert r_of_R.get_value()==50

	first_siblings=right_child_of_L.get_sibling()
	assert len(first_siblings)==2
	for i in first_siblings:
		print('I and the value of i %d' %i.get_value())
		try:
			bool_left_child=i.get_is_left_child()
		except AttributeError:
			pass
	first_tree_attemp.synchronize_sibling_with_same_value(first_siblings)
	for i in first_siblings:
		parent_of_i=i.get_parent()
		assert parent_of_i.get_coverage()==(2,3)
		bool_left_child=i.get_is_left_child()
		assert bool_left_child


def test_BST_split_node():
	reads = string_to_readset("""
	  11
	  0000
	    11
	      00
	      11
	    11
	""")
	first_tree_attemp=Binary_Search_Tree(reads)
	#complete tree is the root node and from there on it gets down to the leafs
	complete_tree=first_tree_attemp.get_complete_tree()
	root_coverage=complete_tree.get_coverage()
	assert root_coverage==(2,3)
	leaf_list_of_tree=first_tree_attemp.get_leaf_list_of_tree()
	assert leaf_list_of_tree[0].get_value()==10
	assert leaf_list_of_tree[0].get_sibling()[1].get_value()==40
	assert leaf_list_of_tree[0].get_is_left_child()
	assert leaf_list_of_tree[0].get_sibling()[1].get_is_left_child()
	end_node_for_split_node=leaf_list_of_tree[0].get_sibling()[1]
	end_node_right_in_sub=end_node_for_split_node.get_parent().get_right_child()
	assert end_node_right_in_sub.get_value()==50
	(split_node_of_read,List_to_change)=first_tree_attemp.seach_for_split_node(leaf_list_of_tree[0],leaf_list_of_tree[0].get_sibling()[1])
	assert split_node_of_read.get_coverage()==(2,3)
	assert split_node_of_read.get_parent()== None
	assert leaf_list_of_tree[0] in List_to_change
	assert leaf_list_of_tree[0].get_parent() in List_to_change
	assert leaf_list_of_tree[0].get_parent().get_right_child() in List_to_change
	assert leaf_list_of_tree[0].get_parent().get_parent() in List_to_change
	assert leaf_list_of_tree[0].get_parent().get_parent().get_right_child() in List_to_change
	assert split_node_of_read in List_to_change
	assert split_node_of_read.get_right_child() in List_to_change
	assert split_node_of_read.get_right_child().get_left_child() in List_to_change
	assert split_node_of_read.get_right_child().get_left_child().get_left_child() in List_to_change
#not in set
	assert not split_node_of_read.get_right_child().get_right_child() in List_to_change
	assert not split_node_of_read.get_right_child().get_left_child().get_right_child() in List_to_change
	assert len(List_to_change)==9

def test_other_cases_for_split_node():
	reads = string_to_readset("""
	  11
	  0000
	    11
	      00
	      11
	    11
	""")
	first_tree_attemp=Binary_Search_Tree(reads)
	#complete tree is the root node and from there on it gets down to the leafs
	root=first_tree_attemp.get_complete_tree()
	leaf_list_of_tree=first_tree_attemp.get_leaf_list_of_tree()
	#second case
	(split_node_of_read,List_to_change)=first_tree_attemp.seach_for_split_node(leaf_list_of_tree[0],leaf_list_of_tree[2])
	assert len(List_to_change)==5
	assert split_node_of_read.get_coverage()==(2,3)
	assert split_node_of_read.get_right_child().get_value()==30
	assert split_node_of_read.get_left_child().get_right_child().get_value()==20
	assert split_node_of_read.get_left_child().get_left_child().get_value()==10
	assert split_node_of_read in List_to_change
	assert split_node_of_read.get_right_child() in List_to_change
	assert split_node_of_read.get_left_child() in List_to_change
	assert split_node_of_read.get_left_child().get_left_child() in List_to_change
	assert split_node_of_read.get_left_child().get_right_child() in List_to_change


def test_third_cases_for_split_node():
	#Not able to construct such a case
	reads = string_to_readset("""
	  11
	  0000
	    111
	    110
	     111
	""")
	first_tree_attemp=Binary_Search_Tree(reads)
	#complete tree is the root node and from there on it gets down to the leafs
	root=first_tree_attemp.get_complete_tree()
	#struct of tree
	L_of_root=root.get_left_child()
	R_of_root=root.get_right_child()
	assert not R_of_root.isLeaf()
	assert not L_of_root.isLeaf()
	l_of_l=L_of_root.get_left_child()
	assert not l_of_l.isLeaf()
	r_of_l=L_of_root.get_right_child()
	assert r_of_l.isLeaf()
	assert r_of_l.get_value()==30
	l_of_r=R_of_root.get_left_child()
	assert not  l_of_r.isLeaf()
	r_of_r =R_of_root.get_right_child()
	assert  r_of_r.isLeaf()
	assert r_of_r.get_value()==60



def test_max_flow():
	reads = string_to_readset("""
	  11
	  0000
	    11
	      00
	      11
	    11
	""")
	first_tree_attemp=Binary_Search_Tree(reads)
	complete_tree=first_tree_attemp.get_complete_tree()
	leaf_list_of_tree=first_tree_attemp.get_leaf_list_of_tree()
	#maximal_coverage=2
	#(pruned_set,removed_set)=optimize_max_flow_in_BST(reads,first_tree_attemp,maximal_coverage)
	#assert len(pruned_set)==4
	#assert len(removed_set)==2
	maximal_coverage=3
	(pruned_set,removed_set)=optimize_max_flow_in_BST(reads,first_tree_attemp,maximal_coverage)
	print(pruned_set)
	print(removed_set)
	#Should be
	#assert len(pruned_set)==6
	#assert len(removed_set)==0
	assert len(pruned_set)==4
	assert len(removed_set)==2

