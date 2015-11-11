import math
from whatshap._core import PyRead, PyReadSet
from whatshap.Binary_Search_Tree import Binary_Search_Tree

# Implementation of the Interval scheduling problem described in the paper of Veli Mäkinen  "Interval scheduling maximizing minimum coverage "


#IMPORTANT TO NOTICE: NOT SUITABLE FOR PAIRED END READS

#First like in the score_based approach and the random approach: Remove the reads which only cover one variant


#TODO:
#Build up a perfect binary search tree with the delimiters of the intervals(so the reads) as leaves.

#Because we need to keep track of the reads which cover the same regions like the reads discarderd (same problem like
# we had in the priority queue) i implemented the Tree structure by myself with not only parent, siblings, which represent the reads covering this node
#Additionally in the initialization each node stores one value its position and later the attributes are added

def reduce_readset_via_max_flow(readset,max_cov):
    tree=Binary_Search_Tree(readset)
    (pruned_set,removed_set)=optimize_max_flow_in_BST(tree,max_cov)
    selected_reads = set()
    #TODO Not sure if the uninformative read count is defined correctly
    #Same as in readselect:
    undecided_reads = set(i for i, read in enumerate(readset) if len(read) >= 2)
    uninformative_read_count=len(readset)-len(undecided_reads)
    for i in pruned_set:
        selected_reads.add(i)
    return selected_reads,uninformative_read_count
    #sliced_reads = reads.subset(selected_reads)





def optimize_max_flow_in_BST(BST, max_cov):
    '''
    Before working with the sibling need to synchronize them and later use the indices of the end nodes for a
    distinct mapping to the reads
    '''
    pruned_for_ending=[]
    removed_for_ending=[]
    leaf_list = BST.get_leaf_list_of_tree()
    for i, leaf in enumerate(leaf_list):
        leaf_value = leaf.get_value()
        siblings = leaf.get_sibling()
        #No need to sort again, done by the construction of leaf list
        only_end_points = [sib for sib in siblings if sib.get_value() > leaf_value]
        BST.synchronize_sibling_with_same_value(only_end_points)
        for end_node in only_end_points:
            (split_node, List_to_change) = BST.seach_for_split_node(leaf, end_node)
            if split_node == None:
                print('Found no split node')
            else:
                split_cov = split_node.get_coverage()
                split_balance = split_node.get_balance()
                selection_criterion = BST.is_crucial(split_cov, max_cov, split_balance)
                if selection_criterion:
                    pruned_for_ending.append(end_node.get_index())
                #Need to change balance fot the involved nodes
                else:
                    removed_for_ending.append(end_node.get_index())
                    step_up_balance(List_to_change)
                    update_till_root(split_node)

    #TODO Need to call a method to select the reads out of the pruned set
    #for test case return both sets

    return (pruned_for_ending,removed_for_ending)

#after each balance step the siblings should be synchronized again, because of all siblings of a node the balance changed
#TODO : Is that needed? Only problem is because in the split node we get all the different siblings of the same value

def step_up_balance(List_to_change):
    #If read is not included need to decrease balance of the nodes connected with them
    for l in List_to_change:
        balance_of_l=l.get_balance()
        l.set_balance(balance_of_l -1)

def update_till_root(split_node):
    '''
    update all coverage till root node depending on the balance
    '''
    #split_node.get_parent()== None means  that split node is the root
    while not(split_node.is_root()):
        split_node_parent=split_node.get_parent()
        min_cov_of_split=split_node.get_min_coverage() +split_node.get_balance()
        max_cov_of_split=split_node.get_max_coverage()+split_node.get_balance()

        if split_node.get_is_left_child():
            other_child=split_node_parent.get_right_child()
            if other_child.isLeaf():
                other_child_coverage_min=other_child.get_coverage()+ other_child.get_balance()
                other_child_coverage_max=other_child.get_coverage()+ other_child.get_balance()
            else:
                other_child_coverage_min=other_child.get_min_coverage()+ other_child.get_balance()
                other_child_coverage_max=other_child.get_max_coverage()+ other_child.get_balance()
        else:
            other_child=split_node_parent.get_left_child()
            if other_child.isLeaf():
                other_child_coverage_min=other_child.get_coverage()+ other_child.get_balance()
                other_child_coverage_max=other_child.get_coverage()+ other_child.get_balance()
            else:
                other_child_coverage_min=other_child.get_min_coverage()+ other_child.get_balance()
                other_child_coverage_max=other_child.get_max_coverage()+ other_child.get_balance()
        parent_split_node_min_coverage=split_node_parent.get_min_coverage()
        parent_split_node_max_coverage=split_node_parent.get_max_coverage()

        if min_cov_of_split<parent_split_node_min_coverage:
            split_node_parent.set_min_coverage(min_cov_of_split)
        if max_cov_of_split>parent_split_node_max_coverage:
            split_node_parent.set_max_coverage(max_cov_of_split)

        split_node=split_node_parent


    #TODO : Look if also include other child coverage_min and max
