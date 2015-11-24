import math
from whatshap._core import PyRead, PyReadSet
from whatshap.Binary_Search_Tree import Binary_Search_Tree
from whatshap.coverage import CovMonitor

# Implementation of the Interval scheduling problem described in the paper of Veli Mäkinen  "Interval scheduling maximizing minimum coverage "


#IMPORTANT TO NOTICE: NOT SUITABLE FOR PAIRED END READS

#First like in the score_based approach and the random approach: Remove the reads which only cover one variant


#TODO:
#Build up a perfect binary search tree with the delimiters of the intervals(so the reads) as leaves.

#Because we need to keep track of the reads which cover the same regions like the reads discarderd (same problem like
# we had in the priority queue) i implemented the Tree structure by myself with not only parent, siblings, which represent the reads covering this node
#Additionally in the initialization each node stores one value its position and later the attributes are added

def look_at_coverage_of_pruned_readset(readset,max_coverage):

    positions = readset.get_positions()
    vcf_indices = {position: index for index, position in enumerate(positions)}
    f = open('Looking_at_coverage', 'w')

    coverages = CovMonitor(len(positions))



    indices_of_reads = set(i for i, read in enumerate(readset) )
    for i in indices_of_reads:
        read_of_i=readset[i]
        print('READ in readset')
        print(read_of_i)
        begin_position=read_of_i[0].position
        begin=vcf_indices.get(begin_position)
        #Not sure if +1 at the end is needed
        #print('read_of_i[len(read_of_i)-1]')
        #print(read_of_i[len(read_of_i)-1])
        end_position=read_of_i[len(read_of_i)-1].position
        end=vcf_indices.get(end_position)
        print("Begin")
        print(begin)
        print("End")
        print(end)
        coverages.add_read(begin, end)
        if coverages.max_coverage_in_range(begin, end) >= max_coverage:
            f.write("coverage exceeded")
            f.write("\t")
            f.write(str(begin_position))
            f.write("\t")
            f.write(str(end_position))
            f.write("\n")
    f.close()


#    for read in readset:
#        print('Read')
#        print(read.getPosition())
#        begin = vcf_indices.get(read.getPosition(0))
#        end = vcf_indices.get(read.getPosition(read.getVariantCount() - 1)) + 1
#        coverages.add_read(begin, end)
#        if coverages.max_coverage_in_range(begin, end) >= max_coverage:
#            f.write("coverage exceeded")
#            f.write("\t")
#            f.write(begin)
#            f.write("\t")
#            f.write(end)
#            f.write("\n")
#    f.close()



def reduce_readset_via_max_flow(readset,max_cov):
    tree=Binary_Search_Tree(readset)
    print("Builded Tree")
    (pruned_set,removed_set)=remove_and_include_reads_from_the_tree(tree,max_cov)
    print('HAS a Pruned and removed set')
    selected_reads = set()
    #TODO Not sure if the uninformative read count is defined correctly
    #Same as in readselect:
    #Already filtered out number of ready by the construction of the tree
    undecided_reads = set(i for i, read in enumerate(readset) if len(read) >= 2)
    uninformative_read_count=len(readset)-len(undecided_reads)
    print('PRUNED SET ')
    print(pruned_set)
    for i in pruned_set:
        #if len(i)>1:
        #    for j in i :
        #       selected_reads.add(j)
        #else:
        #not Neede because i is an integer
        selected_reads.add(i)
    return selected_reads,uninformative_read_count
    #sliced_reads = reads.subset(selected_reads)

def get_all_reads_in_the_range(change_list,start):
    '''
    start= leaf node at this moment
    '''
    print('At calling of get_all_read : Length of change_lsit should be 6 %d' %len(change_list))
    node_combinations=[]
    start_siblings=start.get_sibling()
    already_selected_reads=[index for (sibling,val,index)in start_siblings if val>start.get_value()]
    for i in change_list:
        #if we have a leaf which index is not a delimiter of an already seen read
        if i.isLeaf() :
            print('I value')
            print(i.get_value())
            print(i.get_index())
            if i.get_value()==30:
                print(i)
            for index in i.get_index():
                if index not in already_selected_reads:
                    node_combinations.append((i,index))
    return node_combinations

def remove_and_include_reads_from_the_tree(BST,max_coverage):
    selected_reads=set()
    removed_reads=set()
    crucial_indices=set()
    not_crucial=set()
    leaf_list=BST.get_leaf_list_of_tree()
    #going over the leafs in the tree
    for i,leaf in enumerate(leaf_list):
        print("Looking at Leaf %d" %leaf.get_value())
        leaf_value=leaf.get_value()
        siblings=leaf.get_sibling()
        only_end_points=[[sib,i] for (sib,val,i) in siblings if val> leaf_value]
        #store for every  leaf a list of siblings which xould be removed if coverage is exceeded
        already_seen_siblings=[]
        for end_node in only_end_points:
            #compute split node for every range
            print("Corresponding leaf_node %d" %end_node[0].get_value())
            (split_node,List_to_change,coverage_in_range)=BST.seach_for_split_node(leaf,end_node[0])
            already_seen_siblings.append((split_node,List_to_change,end_node[1]))
            if split_node==None:
                print('Found no split node')
            else:
                split_balance = split_node.get_balance()
                selection_criterion = BST.is_crucial(coverage_in_range, max_coverage, split_balance)
                #This means all reads which cover the read and the end node are crucial, therefore all have to be included
                #in the reduced readset
                if selection_criterion:
                    print('In SELECTION CRITERION  ')
                    print(end_node[1])
                    crucial_indices.add(end_node[1])
                    already_seen_siblings.remove((split_node,List_to_change,end_node[1]))
                    selected_reads.add(end_node[1])
                else:
                    print('NOT crucial')
                    not_crucial.add(end_node[1])


        new_leaf_coverage=leaf.get_coverage() +leaf.get_balance()
        number_need_to_remove=new_leaf_coverage - max_coverage
        print('New leaf coverage %d ' %new_leaf_coverage)
        print('Number of need to remove %d' %number_need_to_remove)
        print("Length of the already selected set %d" %len(already_seen_siblings))
        print('Two sets of pruned and removed and crucial and not crucial ')
        print(selected_reads)
        print(removed_reads)
        print(crucial_indices)
        print(not_crucial)
        #If we have the case that the later reads are crucial but the former were expandable we have to look again at all reads
        # which also cover this region and are expandable


        #First case could delete the number needed from the reads starting at the leaf itself
        if len(already_seen_siblings)>number_need_to_remove:
            print('In first case : already seen big enough')
            while number_need_to_remove >0:
                split,change_list,connecting_index=already_seen_siblings.pop()
                print('Appended to removed_reads %d' %connecting_index)
                removed_reads.add(connecting_index)
                step_up_balance(change_list)
                update_till_root(split)
                number_need_to_remove-=1
                #remainig reads are selcted because they do not interfere with the max coverage
            for (split_n,change_list_n,index) in already_seen_siblings:
                print('Appended index at selected %d '%index)
                selected_reads.add(index)
                not_crucial.add(index)

        #second case : Not enough reads which start at the node there to fullfill the number which has to be removed:
        else:
            if number_need_to_remove==0:
               continue
            print('In second case already seen not enough')
            #remove at first as manny nodes as possible
            while len(already_seen_siblings)!=0:
                print('removing of available in seen')
                split,change_list,connecting_index=already_seen_siblings.pop()
                print('Length of the already selected !=0 %d' %connecting_index)
                removed_reads.add(connecting_index)
                step_up_balance(change_list)
                update_till_root(split)
                number_need_to_remove-=1
            # Getting the reads which end at this node
            only_start_points_in_leaf=[[sib,i] for (sib,val,i) in siblings if val< leaf_value]
            for start_node in only_start_points_in_leaf:
                print('Looking at the other siblings of the leaf where the leaf is the end')
                #so if the read is expandable it could be removed
                if start_node[1] in not_crucial:
                    print('START Node in the only start_points %d '%start_node[0].get_value())
                    print('START Node  index  %d '%start_node[1])
                    (split_node,List_to_change,coverage_in_range)=BST.seach_for_split_node(leaf,start_node[0])
                    removed_reads.add(start_node[1])
                    selected_reads.remove(start_node[1])
                    step_up_balance(List_to_change)
                    update_till_root(split_node)
                    number_need_to_remove-=1
                    #Stop when number of removed reads is fullfilled
                if number_need_to_remove==0:
                    break

            #third case : If still more has to be removed as from this leaf could be done, need to check the range and the nodes there:
            if number_need_to_remove!=0:
                print('In third case nee to look in the ranges')
                big_list_of_changes=set()
                all_siblings_of_leaf=[[sib,i] for (sib,val,i) in siblings]
                other_reads_covering_the_reads=[]
                for sib_nodes in all_siblings_of_leaf:
                    (split_node,List_to_change,coverage_in_range)=BST.seach_for_split_node(leaf,sib_nodes[0])
                    split_balance=split_node.get_balance()
                    #Not seen reasd till there could occure therefor need to check first if crucialot not
                    selection_criterion = BST.is_crucial(coverage_in_range, max_coverage, split_balance)
                    #only if the selection criterion is not fullfilled these reads could be used for the removal
                    if not selection_criterion:
                        big_list_of_changes.union(List_to_change,selection_criterion,sib_nodes[0])
                    else:
                        crucial_indices.add(sib_nodes[1])


                for node in big_list_of_changes:
                    if node.isLeaf():
                        node_value=node.get_value()
                        for (other_siblings,val,i) in node.get_sibling():
                            #by less or bigger the leaf node itself is discarded
                            if (node_value>leaf_value and val<leaf_value and i not in crucial_indices):
                                other_reads_covering_the_reads.append((i,node,other_siblings))
                            if (node_value<leaf_value and val>leaf_value and i not in crucial_indices):
                                other_reads_covering_the_reads.append((i,node,other_siblings))
                while number_need_to_remove!=0:
                    (index,node,other_sibling)=other_reads_covering_the_reads.pop()
                    removed_reads.add(index)
                    #Not sure if needed
                    selected_reads.remove(index)
                    #Search again split node of this couple
                    (split_node,List_to_change,coverage_in_range)=BST.seach_for_split_node()
                    selection_criterion = BST.is_crucial(coverage_in_range, max_coverage, split_balance)
                    if selection_criterion:
                        print("THERE IST SOMETHING WRONG")
                    step_up_balance(List_to_change)
                    update_till_root(split_node)
                    number_need_to_remove-=1



        print('Two sets after while loop of pruned and removed')
        print(selected_reads)
        print(removed_reads)
    return (selected_reads,removed_reads)


def remove_and_include_reads_from_the_tree_former_approach(BST,max_coverage):
    selected_reads=[]
    removed_reads=[]
    crucial_indices=[]
    not_crucial=[]
    leaf_list=BST.get_leaf_list_of_tree()
    #going over the leafs in the tree
    for i,leaf in enumerate(leaf_list):
        print("Looking at Leaf %d" %leaf.get_value())
        print('Coverage of the above leaf %d' %leaf.get_coverage())
        leaf_value=leaf.get_value()
        siblings=leaf.get_sibling()
        only_end_points=[[sib,i] for (sib,val,i) in siblings if val> leaf_value]
        print('Only end points')
        print(only_end_points)
        #store for every  leaf a list of siblings which xould be removed if coverage is exceeded
        already_seen_siblings=[]
        for end_node in only_end_points:
            print('END NODE IN For Loop')
            print(end_node)
            (split_node,List_to_change,coverage_in_range)=BST.seach_for_split_node(leaf,end_node[0])
            print('After calling split node of leaf and end node : Length of change_lsit should be 6 %d' %len(List_to_change))
            print('ASSERTION THAT END NODE IS IN CHANGe LIST Directly after calling search for split node')
            print(end_node[0].get_value())
            print(end_node[0] in List_to_change)

            print('New computed coverage in range ')
            print(coverage_in_range)
            #if end_node[1].isList()
            #if len(end_node[1])>1:
            #    connecting_indices= [j for j in end_node[1] if j in leaf.get_index()]
            #    for k in connecting_indices:
            #        already_seen_siblings.append(split_node,List_to_change,k)
            #else:
            #    already_seen_siblings.append(split_node,List_to_change,end_node[1])
            already_seen_siblings.append((split_node,List_to_change,end_node[1]))
            (s_node,list_to_change,end_node_index) =already_seen_siblings.pop()
            print("Assertions that it is good added at the already seen sibling")
            print(split_node.get_coverage())
            print(end_node[0].get_value())
            print(end_node[0].get_index())
            print(end_node[0] in list_to_change)

            already_seen_siblings.append((split_node,List_to_change,end_node[1]))
            #print('Connecting indices')
            #print(connecting_indices)
            #split the nodes up for the indices

            if split_node==None:
                print('Found no split node')
            else:
                #split_cov = split_node.get_coverage()
                split_balance = split_node.get_balance()
                selection_criterion = BST.is_crucial(coverage_in_range, max_coverage, split_balance)
                #This means all reads which cover the read and the end node are crucial, therefore all have to be included
                #in the reduced readset
                if selection_criterion:
                    print('In SELECTION CRITERION THEREFOR REMOVE OF ')
                    print(end_node[1])
                    crucial_indices.append(end_node[1])
                    already_seen_siblings.remove((split_node,List_to_change,end_node[1]))
                    selected_reads.append(end_node[1])
                   # for k in connecting_indices:
                   #     selected_reads.append(k)



        new_leaf_coverage=leaf.get_coverage() +leaf.get_balance()
        number_need_to_remove=new_leaf_coverage - max_coverage
        print('Before number of need to remove sp where the already selected.pop occures')
        print('Number of need to remove %d' %number_need_to_remove)
        print(new_leaf_coverage)
        print(leaf.get_coverage())
        print(leaf.get_balance())
        print(already_seen_siblings)
        print("Length of the already selected set %d" %len(already_seen_siblings))
        print('Two sets of pruned and removed')
        print(selected_reads)
        print(removed_reads)
        #If we have the case that the later reads are crucial but the former were expandable we have to look again at all reads
        # which also cover this region and are expandable
        if len(already_seen_siblings)<number_need_to_remove:
            print('IF length of already selected lower than number need to remove')
            while len(already_seen_siblings)!=0:
                #Do the same as below till the already seen is empty
                split,change_list,connecting_index=already_seen_siblings.pop()
                print('Length of the already selected !=0' %connecting_index)
                removed_reads.append(connecting_index)
                step_up_balance(change_list)
                update_till_root(split)
            print('ASSERTION THAT END NODE IS IN CHANGe LIST')
            print(end_node[0].get_value())
            print(end_node[0] in change_list)
            print("Connecting_index %d" %connecting_index)
            print(leaf.get_value())
            print(split.get_coverage())
            print("Connecting index of this setting ")
            print(leaf.get_index())
            for i in change_list:
                if i.isLeaf():
                    print("Leaf in change list")
                    print(i.get_value())
            indices=get_all_reads_in_the_range(change_list,leaf)
            print("indices")
            print(indices)
        #The case that the length of the already seen_siblings is 0 so no elements there for removing
        if len(already_seen_siblings)==0 and number_need_to_remove!=0:
            #for the actual leaf:
            only_start_points=[[sib,i] for (sib,val,i) in siblings if val< leaf_value]
            #because all of the possible end points are crucial or removed
            for start_node in only_start_points:
                (split_node,List_to_change,coverage_in_range)=BST.seach_for_split_node(leaf,end_node[0])
                split_balance = split_node.get_balance()
                selection_criterion = BST.is_crucial(coverage_in_range, max_coverage, split_balance)
                if not selection_criterion:
                    removed_reads.append(start_node[1])
                    if start_node[i] in selected_reads:
                        selected_reads.remove(start_node[i])
                        step_up_balance(List_to_change)
                        update_till_root(split_node)
                        number_need_to_remove-=1





        while number_need_to_remove >0:
            #print('Already_seen_siblings')
            #print(already_seen_siblings)
            split,change_list,connecting_index=already_seen_siblings.pop()
            print('Appended to removed_reads %d' %connecting_index)
            removed_reads.append(connecting_index)
            step_up_balance(change_list)
            update_till_root(split)
            number_need_to_remove-=1
        print('Two sets before while loop of pruned and removed')
        print(selected_reads)
        print(removed_reads)
        #Other not decided indices belong in the pruned set
        #print("Before for loop again length of already selected %d" %len(already_seen_siblings))
        for (split_n,change_list_n,index) in already_seen_siblings:
            print('Appended index at selected %d '%index)
            selected_reads.append(index)
            not_crucial.append(index)

        print('Two sets after while loop of pruned and removed')
        print(selected_reads)
        print(removed_reads)
    return (selected_reads,removed_reads)



#Former approach, hopefuly not needed again
def optimize_max_flow_in_BST(BST, max_cov):
    '''
    Before working with the sibling need to synchronize them and later use the indices of the end nodes for a
    distinct mapping to the reads
    '''
    pruned_for_ending=[]
    removed_for_ending=[]
    leaf_list = BST.get_leaf_list_of_tree()
    for i, leaf in enumerate(leaf_list):
        print('Working with Leaf : %d' %leaf.get_value())
        leaf_value = leaf.get_value()
        siblings = leaf.get_sibling()
        #No need to sort again, done by the construction of leaf list
        only_end_points = [sib for (sib,val,index) in siblings if val > leaf_value]
        #for every leaf new list for later selection
        already_selected=[]
        for end_node in only_end_points:
            (split_node, List_to_change) = BST.seach_for_split_node(leaf, end_node)
            connecting_indices=[i for i in end_node.get_index() if i in leaf.get_index()]
            print('Connecting indices between leaf and end node')
            print(connecting_indices)
            already_selected.append((split_node,List_to_change,connecting_indices))
            if split_node == None:
                print('Found no split node')
            else:
                split_cov = split_node.get_coverage()
                split_balance = split_node.get_balance()
                selection_criterion = BST.is_crucial(split_cov, max_cov, split_balance)
                print('Included following end delimiter to the alredy_selected %d ' %end_node.get_value() )
                if selection_criterion:
                    pruned_for_ending.append(end_node.get_index())
                    already_selected.remove((split_node,List_to_change,end_node.get_index()))
                #Need to change balance fot the involved nodes

        new_leaf_coverage=leaf.get_coverage() +leaf.get_balance()
        number_need_to_remove=new_leaf_coverage - max_cov
        print('Before number of need to remove sp where the already selected.pop occures')
        print('Number of need to remove %d' %number_need_to_remove)
        print(new_leaf_coverage)
        print(leaf.get_coverage())
        print(leaf.get_balance())
        print(already_selected)
        print("Length of the already selected set %d" %len(already_selected))
        print('Two sets of pruned and removed')
        print(pruned_for_ending)
        print(removed_for_ending)
        while number_need_to_remove >0:
            #for all not crucial intervals
            #split,change_list,connecting_indices=already_selected[len(already_selected)-1]
            split,change_list,connecting_indices=already_selected.pop()
            ######################former
            #TODO Could happen that we e.g. need to remove 2 nodes and we pop from empty list
            index_of_end=connecting_indices.pop()
            #is a list
            #for index_of_node in leaf.get_index():
            #    if index_of_node in connecting_indices:
            #        removed_for_ending.append(index_of_node)
            #        break
            #    else:
            #        print('SHOULD NOT BE ABLE TO OCCURE')

            ################former
            removed_for_ending.append(index_of_end)
            #print('Call step up balance')
            #print("Value of 40 in change list")
            #help_list=[node.get_value() for node in change_list if node.isLeaf()]
            #help_list2=[node.get_coverage()+node.get_balance() for node in change_list if (node.isLeaf() and node.get_value()==40)]
            #new_set=set(help_list)
            #print(40 in help_list)
            #print(help_list2)
            #print('Leaf %d'%leaf.get_value())
            #print('End_node %d'%index_of_end)
            step_up_balance(change_list)
            update_till_root(split)
            number_need_to_remove-=1
            #new_help_list2=[node.get_coverage()+node.get_balance() for node in change_list if (node.isLeaf() and node.get_value()==40)]
            #print("new_help_list2")
            #print(new_help_list2)
        if len(connecting_indices) !=0:
           for i in connecting_indices:
               pruned_for_ending.append(i)

        print('Two sets before while loop of pruned and removed')
        print(pruned_for_ending)
        print(removed_for_ending)
        #Other not decided indices belong in the pruned set
        print("Before for loop again length of already selected %d" %len(already_selected))
        for (split_n,change_list_n,index) in already_selected:
            pruned_for_ending.append(index.pop())

        print('Two sets after while loop of pruned and removed')
        print(pruned_for_ending)
        print(removed_for_ending)
                #else:
                #    removed_for_ending.append(end_node.get_index())
                #    step_up_balance(List_to_change)
                #    update_till_root(split_node)

    #TODO Need to call a method to select the reads out of the pruned set
    #for test case return both sets

    return (pruned_for_ending,removed_for_ending)

#after each balance step the siblings should be synchronized again, because of all siblings of a node the balance changed
#TODO : Is that needed? Only problem is because in the split node we get all the different siblings of the same value

def step_up_balance(List_to_change):
    print('Step_up_balance')

    #If read is not included need to decrease balance of the nodes connected with them
    for l in List_to_change:
        balance_of_l=l.get_balance()
        l.set_balance(balance_of_l -1)
        #if (l.isLeaf()) and (l.get_value()==50):
        #   print("Leaf with value 50 is found in step up balance")

def update_till_root(split_node):
    '''
    update all coverage till root node depending on the balance
    '''
    print('Update till root')
    #split_node.get_parent()== None means  that split node is the root
    print()
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
