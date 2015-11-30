import math
from whatshap._core import PyRead, PyReadSet
from whatshap.Binary_Search_Tree import Binary_Search_Tree
from whatshap.coverage import CovMonitor


def look_at_coverage_of_pruned_readset(readset,max_coverage):
    '''
    Input: original PyReadSet and the maximal coverage as restricting parameter
    Output: Looks at the coverage at every position in the selected readset to look if the coverage exceeds in the selected readset
    In the reagion where the coverage is exceeded it is written in the file named "Looking_at_coverage"
    '''
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

def reduce_readset_via_max_flow(readset,max_cov):
    '''
    :param: Whole readset  which has to be pruned and the max_cov value of the maximal coverage value
    :return two sets, first the set of the indices which are selected and second the set of the read indices which are not selected
    Starting with building the tree, then decide for every read is it is crucial and then reduce the readset by removing non_crucial
    reads out of the possible set  or adding reads which do not exceed the coverage
    '''
    #Setting up the tree - comments see Binary_Search_Tree class
    tree=Binary_Search_Tree(readset)
    #dividing the reads in crucial and not crucial
    #crucial set is only set of indices, where not_crucial_set is actually a dictionary which stores the index and the  corrsponding
    #delimiters, corresponding split node and the list of nodes, which have to be touched by discarding this node
    (crucial_set,not_crucial_set)=detect_crucial_reads(tree,max_cov)
    print('Crucial set')
    print(crucial_set)
    print('Not crucial set')
    print(not_crucial_set)
    #Should reduce the readset based on the crucial set and the non crucial reads
    #TODO NOT WORKING METHOD AT THE TIME
    (used_set, not_used_set)=reduce_readset_via_crucial_and_non_crucial_reads(tree,max_cov,crucial_set,not_crucial_set)
    #Former approach to remove reads in order to not exceeding the coverage
    #(pruned_set,removed_set)=remove_and_include_reads_from_the_tree(tree,max_cov)

    #Same as in readselect: Need to discard all reads which have length less than 2 , already done in the tree, but later needed
    #for the statistical output therefore computed again here
    undecided_reads = set(i for i, read in enumerate(readset) if len(read) >= 2)
    uninformative_read_count=len(readset)-len(undecided_reads)
    print('PRUNED SET ')
    print(used_set)
    print('Premoved SET ')
    print(not_used_set)
    return used_set,uninformative_read_count


#FORMER
#Needed for remove_and_include_reads_from_the_tree_former_approach
def get_all_reads_in_the_range(change_list,start):
    '''
    :param : change list so the list of nodes depending on a read have to be updatet if the read is not condsidered in the set
    :param: start is the start node of the read
    :return the idea  is to find all leafs in the change_list, so all leaf nodes lying between the start and the end node of a read
    #and only add the indices to the output list, if they are not end nodes, of the current leaf, because these are formerly
    considered in the remove and include reads from the tree method
    start= leaf node at this moment
    '''
    node_combinations=[]
    start_siblings=start.get_sibling()
    already_selected_reads=[index for (sibling,val,index)in start_siblings if val>start.get_value()]
    for i in change_list:
        #if we have a leaf which index is not a delimiter of an already seen read
        if i.isLeaf() :
            for index in i.get_index():
                if index not in already_selected_reads:
                    node_combinations.append((i,index))
    return node_combinations
#FORMER
def remove_and_include_reads_from_the_tree_not_working_approach_with_crucial_set(BST,max_coverage):
    '''
    :param: Getting Tree and the maximal coverage
    :return: (selected_reads,removed_reads)
    The idea is that we have 4 sets, 2 for selected or not selected and 2 for crucial and not crucial
    First thing is to go over the leafs and looking only at the read, where the leaf is the start point.
    Then compute, split node, and looking if the read is crucial or not, if it is crucial the read is in selected read,
    if not, the read is added to a not_crucial set , but also to an already_seen_siblings
    After that looking at the coverage and depending of that if the coverage has to be removed.
    Then need to split between 3 possible scenario :
    1. Have enough reads which start at the leaf which could be removed as the number which has to be removed
        Add remaining to the selected ones
    2. Looking if we do not to reduce the coverage of the node at all
        If only have to reduce it and the already_seen is too small just reduce by the number of elements in the already_seen
        After that look at the reads which end at the leaf and which are already selected, and remove them if they are not crucial
    3. If still the number of reads  which are covering the leaf and need to remove
          Looking at all nodes which are not_crucial and then looking if the node is in selected reads, and if the read in the not crucial is covering
          the leaf and then the read is removed.
    '''
    selected_reads=set()
    removed_reads=set()
    crucial_indices=set()
    not_crucial=set()
    leaf_list=BST.get_leaf_list_of_tree()
    #going over the leafs in the tree
    for i,leaf in enumerate(leaf_list):
        #print("Looking at Leaf %d" %leaf.get_value())
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
                    crucial_indices.add(end_node[1])
                    already_seen_siblings.remove((split_node,List_to_change,end_node[1]))
                    selected_reads.add(end_node[1])
                else:
                    not_crucial.add((end_node[1],end_node[0]))
        new_leaf_coverage=leaf.get_coverage() +leaf.get_balance()
        number_need_to_remove=new_leaf_coverage - max_coverage
        #print('Number of need to remove %d' %number_need_to_remove)
        #print("Length of the already selected set %d" %len(already_seen_siblings))
        #print('Two sets of pruned and removed and crucial and not crucial ')
        #print(selected_reads)
        #print(removed_reads)
        #print(crucial_indices)
        #print(not_crucial)

        #If we have the case that the later reads are crucial but the former were expandable we have to look again at all reads
        # which also cover this region and are expandable

        #First case could delete the number needed from the reads starting at the leaf itself
        if len(already_seen_siblings)>=number_need_to_remove:
            #print('In first case : already seen big enough')
            while number_need_to_remove >0:
                split,change_list,connecting_index=already_seen_siblings.pop()
                #print('Appended to removed_reads %d' %connecting_index)
                removed_reads.add(connecting_index)
                step_up_balance(change_list)
                update_till_root(split)
                number_need_to_remove-=1
                #remainig reads are selcted because they do not interfere with the max coverage
            for (split_n,change_list_n,index) in already_seen_siblings:
                #print('Appended index at selected %d '%index)
                selected_reads.add(index)
                not_crucial.add((index,leaf))

        #second case : Not enough reads which start at the node there to fullfill the number which has to be removed:
        else:
            if number_need_to_remove==0:
               continue
            print('In second case already seen not enough')
            #remove at first as manny nodes as possible
            while len(already_seen_siblings)!=0:
                #print('removing of available in seen')
                split,change_list,connecting_index=already_seen_siblings.pop()
                #print('Length of the already selected !=0 %d' %connecting_index)
                removed_reads.add(connecting_index)
                step_up_balance(change_list)
                update_till_root(split)
                number_need_to_remove-=1
            # Getting the reads which end at this node
            only_start_points_in_leaf=[[sib,i] for (sib,val,i) in siblings if val< leaf_value]
            for start_node in only_start_points_in_leaf:
                #print('Looking at the other siblings of the leaf where the leaf is the end')
                #so if the read is expandable it could be removed
                for (index_not_cruc,node) in not_crucial:
                    if start_node[1] in not_crucial:
                        #print('START Node in the only start_points %d '%start_node[0].get_value())
                        #print('START Node  index  %d '%start_node[1])
                        (split_node,List_to_change,coverage_in_range)=BST.seach_for_split_node(leaf,start_node[0])
                        removed_reads.add(start_node[1])
                        if start_node[1] in selected_reads:
                            #could be somewhere else
                            selected_reads.remove(start_node[1])
                        step_up_balance(List_to_change)
                        update_till_root(split_node)
                        number_need_to_remove-=1
                        #Stop when number of removed reads is fullfilled
                    if number_need_to_remove==0:
                        break

            #third case : If still more has to be removed as from this leaf could be done, need to check the range and the nodes there:

            #other approach for the thrid case :
            #Looking at the not_crucial reads
            while number_need_to_remove!=0:
                checking_set=set()
                for (index,node) in not_crucial:
                    if index in selected_reads:
                        #value_of_not_crucial=node.get_value()
                        node_sibs=node.get_sibling()
                        for (sib,val,i) in node_sibs:
                            if i ==index and val>=leaf_value and node.get_value()<=leaf_value:
                                #print('Covering there with index')
                                (split_node,List_to_change,coverage_in_range)=BST.seach_for_split_node(node,sib)
                                step_up_balance(List_to_change)
                                update_till_root(split_node)
                                number_need_to_remove-=1
                                selected_reads.remove(i)
                                removed_reads.add(i)
                            if i ==index and val<=leaf_value and node.get_value()>=leaf_value:
                                #print('Covering this index')
                                (split_node,List_to_change,coverage_in_range)=BST.seach_for_split_node(node,sib)
                                step_up_balance(List_to_change)
                                update_till_root(split_node)
                                number_need_to_remove-=1
                                selected_reads.remove(i)
                                removed_reads.add(i)
                            if index == i:
                                checking_set.add(index)
                                break

        #print('Two sets after while loop of pruned and removed')
        #print(selected_reads)
        #print(removed_reads)
    return (selected_reads,removed_reads)

def detect_crucial_reads(BST,max_coverage):
    '''
    :param: Getting Tree and max coverage
    :return:2 sets one containing the crucial and one the not_crucial : crucial is set of indices, where not_crucial is a
    dictionary where the keys are the indices and the value the read start and end node , the split node and a list of nodes which
    have to be updated if the read is removed from the original readset
    '''
    crucial_reads=set()
    not_crucial_reads={}
    leaf_list=BST.get_leaf_list_of_tree()
    for leaf in leaf_list:
        siblings_of_leaf=leaf.get_sibling()
        read_startpoints= [(sib,val,i) for (sib,val,i) in siblings_of_leaf  if val> leaf.get_value()]
        for node in read_startpoints:
            (split_node,List_to_change,coverage_in_range)=BST.seach_for_split_node(leaf,node[0])
            split_balance = split_node.get_balance()
            selection_criterion = BST.is_crucial(coverage_in_range, max_coverage, split_balance)
            if selection_criterion:
                #more than just one index ( could occure by double occuring reads, node is the same but index is a list of indices
                #length of node is defined in a leaf as length of the indices of this leaf


                #if len(node[0])>1:
                #    print('In IF and for len of the node ')
                #    print(len(node[0]))
                #    print(node[0])
                #    print(node[1])
                #    print(node[0].get_index())
                #    for i in node[1]:
                #        print(i)
                #        print('CRUCIAL READ')
                #       crucial_reads.add(i)
                #else:
                crucial_reads.add(node[2])
            else:
                not_crucial_reads[node[2]]=((leaf,node,split_node,List_to_change))
    #crucial is a set of leaf start and end point, correspoinding split node and list of nodes, which have to be changed
    #not _crucial is a dictionary, where for every index the same is stored lik in the crucial set
    return (crucial_reads,not_crucial_reads)

#TODO NOt WORKING
def reduce_readset_via_crucial_and_non_crucial_reads(BST, max_coverage,crucial_set,not_crucial_set):
    '''
    :param: Tree, max_coverage, crucial_set of indices and not_crucial_set, corresponding to reads and  split node
    _:return: First set all former selected crucial indices into the selected readset
    Then go over the leaf list : look if we need to remove something :
    Go over the siblings and look if the index is already in the selected set
    '''
    #print('Length of crucial set')
    #print(len(crucial_set))
    #print('Length of not crucial set ')
    #print(len(not_crucial_set))
    selected_set=set()
    not_selected_set=set()
    leaf_list=BST.get_leaf_list_of_tree()
    #first add all indices to the selected which are already crucial
    for i in crucial_set:
        selected_set.add(i)
    #Go over all reads
    for leaf in leaf_list:
        print('At the following leaf %d' %leaf.get_value())
        new_leaf_coverage=leaf.get_coverage() +leaf.get_balance()
        number_need_to_remove=new_leaf_coverage - max_coverage
        print('Number need to remove %d' %number_need_to_remove)
        l_siblings=leaf.get_sibling()
        #get again the reads which start node is the actual leaf
        #while number_need_to_remove>0:
        #    start_nodes= [(sib,i)  for (sib,val,i) in l_siblings if val >leaf.get_value() and i not in selected_set]
        #    for start,index in start_nodes:
        #        print('Number need to remove in for loop %d' %number_need_to_remove)
        #        print('Start nodes')
        #        print('Index')
        #        print(start)
        #        print(index)
        #        #dictionary find results of this index
        #        if index not in selected_set:
        #            ((leaf,node,split_node,List_to_change))=not_crucial_set[index]
        #            print('Not selected the following read')
        #            print(index)
        #            not_selected_set.add(index)
        #            step_up_balance(List_to_change)
        #            update_till_root(split_node)
        #            number_need_to_remove-=1
        #            if number_need_to_remove==0:
        #                print('Found that number need to remove is 0 ')
        #                continue
        #Change the order:
        start_nodes= [(sib,i)  for (sib,val,i) in l_siblings if val >leaf.get_value() and i not in selected_set]
        print('Start nodes list')
        print(start_nodes)
        just_in_this_case_renomved=[]
        #for i in range(0,len(start_nodes)):
        #    (start,index)=start_nodes.pop()
        while number_need_to_remove>0 and len(start_nodes)!=0:
            (start,index)=start_nodes.pop()
            print('Number need to remove in for loop %d' %number_need_to_remove)
            print('Start nodes')
            print('Index')
            print(start)
            print(index)
            #dictionary find results of this index
            if index not in selected_set:
                ((leaf,node,split_node,List_to_change))=not_crucial_set[index]
                print('Not selected the following read')
                print(index)
                not_selected_set.add(index)
                just_in_this_case_renomved.append(index)
                step_up_balance(List_to_change)
                update_till_root(split_node)
                number_need_to_remove-=1
                if number_need_to_remove==0:
                    print('Found that number need to remove is 0 ')

        #for start,index in start_nodes:
        #    print('In for liip over start , index and start_nodes')
        #    while number_need_to_remove>0:
                #print('Number need to remove in for loop %d' %number_need_to_remove)
                #print('Start nodes')
                #print('Index')
                #print(start)
                #print(index)
                #dictionary find results of this index
                #if index not in selected_set:
                #    ((leaf,node,split_node,List_to_change))=not_crucial_set[index]
                #    print('Not selected the following read')
                #    print(index)
                #    not_selected_set.add(index)
                #    just_in_this_case_renomved.append(index)
                #    step_up_balance(List_to_change)
                #    update_till_root(split_node)
                 #   number_need_to_remove-=1
                 #   if number_need_to_remove==0:
                 #       print('Found that number need to remove is 0 ')
        print('Just in this case removed')
        print(just_in_this_case_renomved)
        if start_nodes==[] and number_need_to_remove>0:
            print('No start nodes available but need to remove more')
            print('Not crucial keys')
            print(not_crucial_set.keys())
            not_crucial_indices=set(not_crucial_set.keys())
            print(not_crucial_indices)
            intersection_between_not_crucial_and_already_selected=not_crucial_indices.intersection(selected_set)
            print('intersection_between_not_crucial_and_already_selected')
            print(intersection_between_not_crucial_and_already_selected)
            #Go over this intersecdtion set and look which cover the leaf
            # while number_need_to_remove>0:
            for index in intersection_between_not_crucial_and_already_selected:
                ((start_node,end_node,split_node,List_to_change))=not_crucial_set[index]
                #if actual leaf is influenced by this read
                if leaf in List_to_change and number_need_to_remove>0 and index not in just_in_this_case_renomved:
                    print('Index  where leaf is in List to change %d' %index)
                    not_selected_set.add(index)
                    selected_set.remove(index)
                    step_up_balance(List_to_change)
                    update_till_root(split_node)
                    number_need_to_remove-=1





        #TODO Maybe endless loop
        if number_need_to_remove==0:
            adding_indices=[i for (sib,val,i) in l_siblings if val > leaf.get_value() and i not in not_selected_set]
            #In this loop
            print(' In Adding indices here')
            for index in adding_indices:
                print('Index')
                print(index)
                selected_set.add(index)
        print('TWO sets selected and not selected')
        print(selected_set)
        print(not_selected_set)
    return (selected_set,not_selected_set)



#Former Approach with slightly differnt third step
def remove_and_include_reads_from_the_tree(BST,max_coverage):
    '''
    :param: Tree and maximal coverage
    :return: As the other method before returning selected reads and removed reads
    Till end of step 2 is  the same as inr remove_and_include_reads_from_the_tree_ not_working_approach_with_crucial_set
    Step 3:
    Go over all siblings of the leaf and union the corresponding List_to_change which is computed by the split node search
    This should help, by detecing all leaf nodes, which should be all nodes, which lead to a major coverage of the leaf
    so for every leaf node excluding the leaf and its siblings, and detect reads which lie over the leaf
    Reduce the readset by these nodes .
    '''
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
        if len(already_seen_siblings)>=number_need_to_remove:
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
                    if start_node[1] in selected_reads:
                        #could be somewhere else
                        selected_reads.remove(start_node[1])
                    step_up_balance(List_to_change)
                    update_till_root(split_node)
                    number_need_to_remove-=1
                    #Stop when number of removed reads is fullfilled
                if number_need_to_remove==0:
                    break

            #third case : If still more has to be removed as from this leaf could be done, need to check the range and the nodes there:
            if number_need_to_remove!=0:
                print('In third case need to look in the ranges')
                big_list_of_changes=set()
                all_siblings_of_leaf=[[sib,i] for (sib,val,i) in siblings]
                other_reads_covering_the_reads=[]
                for sib_nodes in all_siblings_of_leaf:
                    print("sib_nodes for the finding of other siblings of the leaf nodes in this range")
                    print(sib_nodes[0].get_value())
                    print(sib_nodes[1])
                    (split_node,List_to_change,coverage_in_range)=BST.seach_for_split_node(leaf,sib_nodes[0])
                    split_balance=split_node.get_balance()
                    #Not seen reasd till there could occure therefor need to check first if crucialot not
                    selection_criterion = BST.is_crucial(coverage_in_range, max_coverage, split_balance)
                    #only if the selection criterion is not fullfilled these reads could be used for the removal
                    if not selection_criterion:
                        print('In selection criterion')
                        print('List_to_change which is there ')
                        print(List_to_change)
                        big_list_of_changes=big_list_of_changes.union(List_to_change)
                        print(len(big_list_of_changes))
                    else:
                        print('Not crucial')
                        crucial_indices.add(sib_nodes[1])

                print("len(big_list_of_changes)")
                print(len(big_list_of_changes))
                for node in big_list_of_changes:
                    print('In For loop')
                    if node.isLeaf():
                        print("NODES HAVE TO BE LEAF ")
                        node_value=node.get_value()
                        for (other_siblings,val,i) in node.get_sibling():
                            #by less or bigger the leaf node itself is discarded
                            if (node_value>leaf_value and val<leaf_value and i not in crucial_indices):
                                print('Added in other_reads_ in_first if ')
                                other_reads_covering_the_reads.append((i,node,other_siblings))
                            if (node_value<leaf_value and val>leaf_value and i not in crucial_indices):
                                print('Added in other_reads_ in_second if ')
                                other_reads_covering_the_reads.append((i,node,other_siblings))
                if len(other_reads_covering_the_reads)==0:
                    print("Nothing to remove")
                    print("Number needs to remove %d" %number_need_to_remove)
                while number_need_to_remove!=0:
                    (index,node,other_sibling)=other_reads_covering_the_reads.pop()
                    removed_reads.add(index)
                    #Not sure if needed
                    if index in selected_reads:
                        selected_reads.remove(index)
                    #Search again split node of this couple
                    (split_node,List_to_change,coverage_in_range)=BST.seach_for_split_node(node,other_sibling)
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
    '''
    :param:BST and max_coverage
    :return Like in the other methods starting with going over the leafs
    Getting the nodes which are start at the leafs. Store the found node, with split node and change list in a list
    And then for every element of this list  assert if the read is crucial or not.



    '''
    selected_reads=[]
    removed_reads=[]
    crucial_indices=[]
    not_crucial=[]
    leaf_list=BST.get_leaf_list_of_tree()
    #going over the leafs in the tree
    for i,leaf in enumerate(leaf_list):
        #print("Looking at Leaf %d" %leaf.get_value())
        #print('Coverage of the above leaf %d' %leaf.get_coverage())
        leaf_value=leaf.get_value()
        siblings=leaf.get_sibling()
        only_end_points=[[sib,i] for (sib,val,i) in siblings if val> leaf_value]
        #print('Only end points')
        #print(only_end_points)
        #store for every  leaf a list of siblings which xould be removed if coverage is exceeded
        already_seen_siblings=[]
        for end_node in only_end_points:
            #print('END NODE IN For Loop')
            #print(end_node)
            (split_node,List_to_change,coverage_in_range)=BST.seach_for_split_node(leaf,end_node[0])
            #print('After calling split node of leaf and end node : Length of change_lsit should be 6 %d' %len(List_to_change))
            #print('ASSERTION THAT END NODE IS IN CHANGe LIST Directly after calling search for split node')
            #print(end_node[0].get_value())
            #print(end_node[0] in List_to_change)

            #print('New computed coverage in range ')
            #print(coverage_in_range)
            already_seen_siblings.append((split_node,List_to_change,end_node[1]))

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
