from phasingutils import string_to_readset
from whatshap.maxflow import  reduce_readset_via_max_flow
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
	only_end_points=[sib for (sib,value,index) in siblings if value>leaf_value]

	assert only_end_points[0].get_value() ==30
	assert only_end_points[1].get_value() ==30
	(split_node_of_read,List_to_change,coverage_range)=tree.seach_for_split_node(leaf_list_of_tree[0],only_end_points[0])
	assert split_node_of_read.get_coverage()==(2,3)
	assert split_node_of_read.get_right_child().get_value()==30
	assert len(List_to_change)==3
	#Geht nicht mehr da List to change nun set ist
	#assert List_to_change[0]==split_node_of_read
	assert split_node_of_read in List_to_change
	(split_node_of_read,List_to_change,coverage_range)=tree.seach_for_split_node(leaf_list_of_tree[0],only_end_points[1])
	assert split_node_of_read.get_coverage()==(2,3)
	assert split_node_of_read.get_right_child().get_value()==30
	leaf_value=leaf_list_of_tree[1].get_value()
	siblings=leaf_list_of_tree[1].get_sibling()
	only_end_points=[sib for (sib,value,index) in siblings if value>leaf_value]
	assert only_end_points[0].get_value() ==40
	(split_node_of_read,List_to_change,coverage_range)=tree.seach_for_split_node(leaf_list_of_tree[1],only_end_points[0])
	assert split_node_of_read.get_coverage()==(2,3)
	assert split_node_of_read==root_node

def test_reveal_whole_tree_structure_for_simple_case():
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
	first_of_value_30=siblings_of_first_leaf[0][0]
	second_of_value_30=siblings_of_first_leaf[1][0]
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
	assert 0 in first_of_value_30.get_index()
	assert 1 in second_of_value_30.get_index()
	#assert get_coverage
	assert first_of_value_30.get_coverage()==3
	assert second_of_value_30.get_coverage()==3
	second_leaf=leaf_list_of_tree[1]
	assert second_leaf.get_value()==30
	assert not second_leaf.get_is_left_child()
	siblings_of_second_leaf=second_leaf.get_sibling()
	assert len(siblings_of_second_leaf)==3
	assert len(siblings_of_second_leaf)==3
	first_of_value=siblings_of_second_leaf[0][0]
	second_of_value=siblings_of_second_leaf[1][0]
	third_of_value=siblings_of_second_leaf[2][0]
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

	assert first_of_value.get_index()==[0,1]
	assert  1 in second_of_value.get_index()
	assert 2 in third_of_value.get_index()
	#assert get_coverage
	assert first_of_value.get_coverage()==2
	assert second_of_value.get_coverage()==2
	assert third_of_value.get_coverage()==3
	third_leaf=leaf_list_of_tree[2]
	assert third_leaf.get_value()==40
	assert third_leaf.get_is_left_child()
	siblings_of_third_leaf=third_leaf.get_sibling()
	assert len(siblings_of_third_leaf)==3
	first_of_value=siblings_of_third_leaf[0][0]
	second_of_value=siblings_of_third_leaf[1][0]
	third_of_value=siblings_of_third_leaf[2][0]
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
	assert 2 in first_of_value.get_index()
	assert 3 in second_of_value.get_index()
	assert 4 in third_of_value.get_index()
	#assert get_coverage
	assert first_of_value.get_coverage()==3
	assert second_of_value.get_coverage()==2
	assert third_of_value.get_coverage()==2
	forth_leaf=leaf_list_of_tree[3]
	assert forth_leaf.get_value()==50
	assert not forth_leaf.get_is_left_child()
	siblings_of_forth_leaf=forth_leaf.get_sibling()
	assert len(siblings_of_forth_leaf)==2
	first_of_value=siblings_of_forth_leaf[0][0]
	second_of_value=siblings_of_forth_leaf[1][0]
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
	assert 3 in first_of_value.get_index()
	assert 4 in second_of_value.get_index()
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
	only_end_points=[sib for (sib,value,index) in siblings if value>leaf_value]
	assert only_end_points[0].get_value() ==30
	assert only_end_points[1].get_value() ==30
	(split_node_of_read,List_to_change,coverage_range)=tree.seach_for_split_node(leaf_list_of_tree[0],only_end_points[0])
	#List to change should include 3 nodes, the 2 nodes and the parent node of both
	assert len(List_to_change) ==3
	out_of_list_to_set=set(List_to_change)
	assert leaf_list_of_tree[0] in out_of_list_to_set
	assert only_end_points[0] in out_of_list_to_set
	assert leaf_list_of_tree[0].get_parent() in out_of_list_to_set

	(split_node_of_read,List_to_change,coverage_range)=tree.seach_for_split_node(leaf_list_of_tree[0],only_end_points[1])
	assert len(List_to_change) ==3
	out_of_list_to_set=set(List_to_change)
	assert leaf_list_of_tree[0] in out_of_list_to_set
	assert only_end_points[1] in out_of_list_to_set
	assert only_end_points[1].get_parent() in out_of_list_to_set
	assert leaf_list_of_tree[0].get_parent() in out_of_list_to_set


	middle_leaf=leaf_list_of_tree[1]
	siblings=middle_leaf.get_sibling()
	assert middle_leaf.get_value()==30
	only_end_points=[sib for (sib,val,index) in siblings if val>leaf_value]
	assert len(only_end_points)==1
	assert only_end_points[0].get_value()==40
	(split_node_of_read,List_to_change,coverage_range)=tree.seach_for_split_node(middle_leaf,only_end_points[0])
	assert split_node_of_read==middle_leaf.get_parent().get_parent()
	print("Length of list to change : %d" %len(List_to_change))
	print(List_to_change)
	assert len(List_to_change)==5
	out_of_list_to_set=set(List_to_change)
	assert middle_leaf in out_of_list_to_set
	assert middle_leaf.get_parent() in out_of_list_to_set
	assert only_end_points[0] in out_of_list_to_set
	assert only_end_points[0].get_parent() in out_of_list_to_set
	assert split_node_of_read in out_of_list_to_set



def test_new_struct_of_siblings():
	reads = string_to_readset("""
	  11
	  0000
	    11
	      00
	      11
	    11
	""")
	first_new_tree=Binary_Search_Tree(reads)
	leaf_after_first_datatset=first_new_tree.get_leaf_list_of_tree()

	node=leaf_after_first_datatset[0]
	assert node.get_coverage()==2
	assert node.get_value() ==10
	assert len(node.get_sibling())==2
	read_index=node.get_sibling()[0][2]
	read_sibling_val=node.get_sibling()[0][1]
	assert read_sibling_val==20
	assert read_index==0
	node=leaf_after_first_datatset[1]
	assert node.get_coverage()==2
	assert node.get_value() ==20
	assert len(node.get_sibling())==1
	read_index=node.get_sibling()[0][2]
	read_sibling_val=node.get_sibling()[0][1]
	assert read_sibling_val==10
	assert read_index==0
	node=leaf_after_first_datatset[2]
	assert node.get_coverage()==3
	assert node.get_value() ==30
	assert len(node.get_sibling())==2
	read_index=node.get_sibling()[0][2]
	read_sibling_val=node.get_sibling()[0][1]
	assert read_sibling_val==40
	assert read_index==2
	node=leaf_after_first_datatset[3]
	assert node.get_coverage()==3
	assert node.get_value() ==40
	assert len(node.get_sibling())==3
	read_index=node.get_sibling()[0][2]
	read_sibling_val=node.get_sibling()[0][1]
	assert read_sibling_val==10
	assert read_index==1
	node=leaf_after_first_datatset[4]
	assert node.get_coverage()==2
	assert node.get_value() ==50
	assert len(node.get_sibling())==2
	read_index=node.get_sibling()[0][2]
	read_sibling_val=node.get_sibling()[0][1]
	assert read_sibling_val==60
	assert read_index==3
	node=leaf_after_first_datatset[5]
	assert node.get_coverage()==2
	assert node.get_value() ==60
	assert len(node.get_sibling())==2
	read_index=node.get_sibling()[0][2]
	read_sibling_val=node.get_sibling()[0][1]
	assert read_sibling_val==50
	assert read_index==3
	read_index=node.get_sibling()[1][2]
	read_sibling_val=node.get_sibling()[1][1]
	assert read_sibling_val==50
	assert read_index==4
	value=first_new_tree.bottom_up_construction_of_BST(leaf_after_first_datatset)
	assert len(value)==4
	inner_node_list=first_new_tree.build_balanced_binary_tree(value,[])
	assert len(inner_node_list)==3


def test_index_of_reads():
	reads = string_to_readset("""
	  11
	  0000
	    11
	      00
	      11
	    11
	""")
	tree=Binary_Search_Tree(reads)
	leaf_list=tree.get_leaf_list_of_tree()
	first_index=leaf_list[0].get_index()
	assert len(first_index)==2
	assert 0 in first_index
	assert 1 in first_index
	second_index=leaf_list[1].get_index()
	assert len(second_index)==1
	assert 0 in second_index
	third_index=leaf_list[2].get_index()
	assert len(third_index)==2
	assert 2 in third_index
	assert 5 in third_index
	forth_index=leaf_list[3].get_index()
	assert len(forth_index)==3
	assert 1 in forth_index
	assert 2 in forth_index
	assert 5 in forth_index
	fifth_index=leaf_list[4].get_index()
	assert len(fifth_index)==2
	assert 3 in fifth_index
	assert 4 in fifth_index
	six_index=leaf_list[5].get_index()
	assert len(six_index)==2
	assert 3 in six_index
	assert 4 in six_index




####Working tests till here
#





def test_simple_case_max_flow():
	reads=string_to_readset("""
 	 111
 	 000
 	   11
 	    00
 	    11
 	""")#
	tree=Binary_Search_Tree(reads)
	root_node=tree.get_complete_tree()
	leaf_list_of_tree=tree.get_leaf_list_of_tree()
	maximal_coverage=2
	pruned_set,indices=reduce_readset_via_max_flow(reads,maximal_coverage)
	#print('Pruned SEt')
 	#print(pruned_set)
 	#print('Removed set')
 	#print(removed_set)
	assert len(pruned_set)==4
	#assert len(removed_set)==1


def test_max_flow_1_unique():
	reads = string_to_readset("""
 	  11
 	  0000
 	    11
 	      00
 	      11
 	    11
 	""")
	first_tree_attemp=Binary_Search_Tree(reads)
	maximal_coverage=1
	(pruned_set,indices)=reduce_readset_via_max_flow(reads,maximal_coverage)
	assert len(pruned_set)==3
	#assert len(removed_set)==3



def test_simple_case_max_flow_2():
	reads=string_to_readset("""
 	 111
 	 000
 	   11
 	    00
 	    11
 	""")
	tree=Binary_Search_Tree(reads)
	maximal_coverage=1
	pruned_set,indices=reduce_readset_via_max_flow(reads,maximal_coverage)
	assert len(pruned_set)==2
	#assert len(removed_set)==3


# 	#Is working

def test_simple_case_max_flow_3():
	reads=string_to_readset("""
	 111
	 000
	   11
	    00
	    11
	""")
	tree=Binary_Search_Tree(reads)
	maximal_coverage=3
	pruned_set,indices=reduce_readset_via_max_flow(reads,maximal_coverage)
	assert len(pruned_set)==5
	#assert len(removed_set)==0

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
	assert left_child.get_coverage()==(1,3)
	assert right_child.get_coverage()==(1,2)
	assert left_child.get_balance()==0
	assert right_child.get_balance()==0
	assert left_child.get_min_coverage()==1
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
	assert R_of_root.get_coverage()==(2,2)
	left_child_of_L=L_of_root.get_left_child()
	right_child_of_L=L_of_root.get_right_child()
	left_child_of_R=R_of_root.get_left_child()
	right_child_of_R=R_of_root.get_right_child()
	assert not right_child_of_L.isLeaf()
	assert not left_child_of_L.isLeaf()
	assert right_child_of_R.isLeaf()
	assert right_child_of_R.get_value()==60
	assert left_child_of_R.isLeaf()
	assert left_child_of_R.get_value()==50
	l_of_L=left_child_of_L.get_left_child()
	r_of_L=left_child_of_L.get_right_child()
	third=right_child_of_L.get_left_child()
	forth=right_child_of_L.get_right_child()
	assert l_of_L.get_value()==10
	assert r_of_L.get_value()==20
	assert third.get_value()==30
	assert forth.get_value()==40
	first_siblings=l_of_L.get_sibling()
	assert len(first_siblings)==2
	(node,value,index)=first_siblings[0]
	assert value==20
	assert index==0
	(node_2,value_2,index_2)=first_siblings[1]
	assert value_2==40
	assert index_2==1
	second_siblings=r_of_L.get_sibling()
	assert len(second_siblings)==1
	(node,value,index)=second_siblings[0]
	assert value==10
	assert index==0
	third_siblings=third.get_sibling()
	assert len(third_siblings)==2
	(node,value,index)=third_siblings[0]
	assert value==40

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
	assert leaf_list_of_tree[0].get_sibling()[1][0].get_value()==40
	assert leaf_list_of_tree[0].get_is_left_child()
	assert leaf_list_of_tree[0].get_sibling()[1][0].get_is_right_child()
	end_node_for_split_node=leaf_list_of_tree[0].get_sibling()[1][0]
	end_node_left_in_sub=end_node_for_split_node.get_parent().get_left_child()
	assert end_node_left_in_sub.get_value()==30
	(split_node_of_read,List_to_change,coverage_range)=first_tree_attemp.seach_for_split_node(leaf_list_of_tree[0],leaf_list_of_tree[0].get_sibling()[1][0])
	assert split_node_of_read.get_coverage()==(2,3)
	assert not split_node_of_read.get_parent()== None
	#root node
	assert split_node_of_read.get_parent().get_parent()== None
	assert leaf_list_of_tree[0] in List_to_change
	assert leaf_list_of_tree[0].get_parent() in List_to_change
	assert leaf_list_of_tree[0].get_parent().get_right_child() in List_to_change
	assert leaf_list_of_tree[0].get_parent().get_parent() in List_to_change
	assert leaf_list_of_tree[0].get_parent().get_parent().get_right_child() in List_to_change
	assert split_node_of_read in List_to_change
	assert split_node_of_read.get_right_child() in List_to_change
	assert split_node_of_read.get_right_child().get_left_child() in List_to_change
#not in set
	assert split_node_of_read.get_right_child().get_right_child() in List_to_change
	assert len(List_to_change)==7


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
	assert leaf_list_of_tree[0].get_value()==10
	assert leaf_list_of_tree[2].get_value()==30
	(split_node_of_read,List_to_change,coverage_range)=first_tree_attemp.seach_for_split_node(leaf_list_of_tree[0],leaf_list_of_tree[2])
	print("Length of list to change : %d" %len(List_to_change))
	assert len(List_to_change)==6
	assert split_node_of_read.get_coverage()==(2,3)
	assert split_node_of_read.get_right_child().get_left_child().get_value()==30
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
	assert not r_of_l.isLeaf()
	l_of_r=R_of_root.get_left_child()
	assert  l_of_r.isLeaf()
	assert  l_of_r.get_value()==50
	r_of_r =R_of_root.get_right_child()
	assert  r_of_r.isLeaf()
	assert r_of_r.get_value()==60

def test_third_case_for_split_node_2():
	reads = string_to_readset("""
	  11
	  0000
	    11
	      00
	      11
	    11
	""")
	first_tree_attemp=Binary_Search_Tree(reads)
	leaf_list_of_tree=first_tree_attemp.get_leaf_list_of_tree()
	(split_node,List_to_change,coverage_range)=first_tree_attemp.seach_for_split_node(leaf_list_of_tree[2],leaf_list_of_tree[3])
	assert len(List_to_change)==3
	assert not leaf_list_of_tree[4] in List_to_change
	assert leaf_list_of_tree[2] in List_to_change
	assert leaf_list_of_tree[3] in List_to_change
	root=first_tree_attemp.get_complete_tree()
	assert not root in List_to_change
	assert not root.get_right_child() in List_to_change
	assert not root.get_left_child() in List_to_change
	assert root.get_left_child().get_right_child() in List_to_change


def test_max_flow_1():
	reads = string_to_readset("""
	  11
	  0000
	    11
	      00
	      11
	    11
	""")
	first_tree_attemp=Binary_Search_Tree(reads)
	maximal_coverage=1
	(pruned_set,indices)=reduce_readset_via_max_flow(reads,maximal_coverage)
	assert len(pruned_set)==3
	#assert len(removed_set)==3



def test_max_flow_2():
	reads = string_to_readset("""
	  11
	  0000
	    11
	      00
	      11
	    11
	""")
	first_tree_attemp=Binary_Search_Tree(reads)
	maximal_coverage=2
	(pruned_set,indices)=reduce_readset_via_max_flow(reads,maximal_coverage)
	assert len(pruned_set)==5
	#assert len(removed_set)==1



def test_max_flow_3():
	reads = string_to_readset("""
	  11
	  0000
	    11
	      00
	      11
	    11
	""")
	first_tree_attemp=Binary_Search_Tree(reads)
	maximal_coverage=3
	(pruned_set,removed_set)=reduce_readset_via_max_flow(reads,maximal_coverage)
	print('pruned_set')
	print(pruned_set)
	assert len(pruned_set)==6
	#assert len(removed_set)==0



def test_split_node_second_case():
	reads = string_to_readset("""
	  11
	  0000
	    11
	      00
	      11
	    11
	""")
	tree=Binary_Search_Tree(reads)
	leaf_list_of_tree=tree.get_leaf_list_of_tree()
	#go over leafs =reads

	leaf_value=leaf_list_of_tree[0].get_value()
	siblings=leaf_list_of_tree[0].get_sibling()
	only_end_points=[sib for (sib,value,index) in siblings if value>leaf_value]

	assert leaf_list_of_tree[3].get_value() ==40
	assert leaf_list_of_tree[4].get_value() ==50
	(split_node_of_read,List_to_change,coverage_range)=tree.seach_for_split_node(leaf_list_of_tree[3],leaf_list_of_tree[4])
	print(List_to_change)
	print(len(List_to_change))
	assert len(List_to_change)==6
	assert split_node_of_read in List_to_change
	assert leaf_list_of_tree[3] in List_to_change
	assert leaf_list_of_tree[4] in List_to_change
	assert leaf_list_of_tree[4].get_parent() in List_to_change
	assert leaf_list_of_tree[3].get_parent() in List_to_change
	assert leaf_list_of_tree[3].get_parent().get_parent() in List_to_change
	print(leaf_list_of_tree[3].get_parent().get_left_child().get_value())
	assert not leaf_list_of_tree[3].get_parent().get_left_child() in List_to_change
	assert not leaf_list_of_tree[4].get_parent().get_right_child() in List_to_change
	assert split_node_of_read==tree.get_complete_tree()


def test_whole_algorithm ():
	reads = string_to_readset("""
	  11
	  0000
	    11
	      00
	      11
	    11
	""")
	maximum_coverage=2
	selected_reads,uninformative_read_count=reduce_readset_via_max_flow(reads,maximum_coverage)
	assert len(selected_reads)==5



def test_is_crucial_of_maxflow():
	reads = string_to_readset("""
	  11
	  0000
	    111
	    110
	     111
	""")
	maximum_coverage=2
	tree=Binary_Search_Tree(reads)
	leaf_list=tree.get_leaf_list_of_tree()
	assert leaf_list[2].get_value()==30
	assert leaf_list[4].get_value()==50
	split_node,List_of_nodes,coverage_of_range=tree.seach_for_split_node(leaf_list[2],leaf_list[4])
	assert coverage_of_range==(3,4)
	assert split_node.get_coverage()==(1,4)

	split_node,List_of_nodes,coverage_of_range=tree.seach_for_split_node(leaf_list[0],leaf_list[4])
	assert coverage_of_range==(2,4)
	assert split_node.get_coverage()==(1,4)

	split_node,List_of_nodes,coverage_of_range=tree.seach_for_split_node(leaf_list[1],leaf_list[2])
	assert coverage_of_range==(2,3)
	assert split_node.get_coverage()==(2,4)

	split_node,List_of_nodes,coverage_of_range=tree.seach_for_split_node(leaf_list[3],leaf_list[5])
	assert coverage_of_range==(1,4)
	assert split_node.get_coverage()==(1,4)

	assert leaf_list[3].get_parent() in List_of_nodes
	assert leaf_list[3].get_parent().get_parent() in List_of_nodes
	assert not leaf_list[3].get_parent().get_left_child() in List_of_nodes
	assert leaf_list[5].get_parent() in List_of_nodes
	assert leaf_list[5] in List_of_nodes
	assert leaf_list[5].get_parent().get_left_child()in List_of_nodes




#TODO NOt working because wrong assumption in the  paper
def test_whole_algorithm_for_simple_case():#
	#Not able to construct such a case
	reads = string_to_readset("""
	  11
	  0000
	    111
	    110
	     111
	""")
	maximum_coverage=2
	selected_reads,uninformative_read_count=reduce_readset_via_max_flow(reads,maximum_coverage)
	assert len(selected_reads)==3

def test_split_node_in_this_case():#
	#Not able to construct such a case
	reads = string_to_readset("""
	  11
	  0000
	    111
	    110
	     111
	""")
	maximum_coverage=2
	tree=Binary_Search_Tree(reads)
	leaf_list=tree.get_leaf_list_of_tree()
	split_node,List_of_nodes,coverage_of_range=tree.seach_for_split_node(leaf_list[3],leaf_list[5])
	assert coverage_of_range==(1,4)
	assert split_node.get_coverage()==(1,4)
	assert leaf_list[3].get_parent() in List_of_nodes
	assert leaf_list[3].get_parent().get_parent() in List_of_nodes
	assert not leaf_list[3].get_parent().get_left_child() in List_of_nodes
	assert leaf_list[5].get_parent() in List_of_nodes
	assert leaf_list[5] in List_of_nodes
	assert leaf_list[5].get_parent().get_left_child()in List_of_nodes
	assert len(List_of_nodes)==7

def test_third_case_of_the_selection_of_reads():
	reads = string_to_readset("""
	  111
	  000
	    11
	  00000
	  11111
		""")
	maximum_coverage=3
	selected_reads,uninformative_read_count=reduce_readset_via_max_flow(reads,maximum_coverage)
	assert len(selected_reads)==3


def test_case_for_using_range():
	reads=string_to_readset("""
	111111
	000000
	111111
	    11
	    000
	""")
	maximum_coverage=3
	selected_reads,uninformative_read_count=reduce_readset_via_max_flow(reads,maximum_coverage)
	assert len(selected_reads)==3
