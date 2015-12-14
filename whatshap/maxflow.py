import math
from whatshap._core import PyRead, PyReadSet
from whatshap.Binary_Search_Tree import Binary_Search_Tree
from whatshap.coverage import CovMonitor

#Degbugging
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
        #print('READ in readset')
        #print(read_of_i)
        begin_position=read_of_i[0].position
        begin=vcf_indices.get(begin_position)
        end_position=read_of_i[len(read_of_i)-1].position
        end=vcf_indices.get(end_position)
        #print("Begin")
        #print(begin)
        #print("End")
        #print(end)
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


    undecided_reads = set(i for i, read in enumerate(readset) if len(read) >= 2)
    #Assertion that all reads are assigned to either one of the sets/dictionary
    assert (len(crucial_set) + len(not_crucial_set))==len(undecided_reads)
    not_crucial_keys=set(not_crucial_set.keys())
    intersect_of_crucial_and_not_crucial=not_crucial_keys.intersection(crucial_set)
    #Assertion that no index has a double occurence
    assert len(intersect_of_crucial_and_not_crucial)==0
    #explore_crucial_set(extended_crucial_set,tree,max_cov)

    #Should reduce the readset based on the crucial set and the non crucial reads
    (used_set, not_used_set)=reduce_readset_via_crucial_and_non_crucial_reads(tree,max_cov,crucial_set,not_crucial_set)
    #Former approach to remove reads in order to not exceeding the coverage
    #(pruned_set,removed_set)=remove_and_include_reads_from_the_tree(tree,max_cov)

    #Same as in readselect: Need to discard all reads which have length less than 2 , already done in the tree, but later needed
    #for the statistical output therefore computed again here
    uninformative_read_count=len(readset)-len(undecided_reads)
    return used_set,uninformative_read_count

#only for debugging
def explore_crucial_set(crucial_set,tree,max_cov):
    leaf_list=tree.get_leaf_list_of_tree()
    for (start_node,end_node,split,List_to_change,index) in crucial_set.values():
        s_val=start_node.get_value()
        e_val=end_node[0].get_value()
        if start_node.isLeaf():
            print('Start node is leaf')
        else:
            print('Start node is not leaf')
        if end_node[0].isLeaf():
            print('End node is leaf')
        else:
            print('End node is not leaf')
        if index==23710:
            if start_node.isLeaf():
                print('start_node.get_value()')
                print(start_node.get_value())
            if end_node[0].isLeaf():
                print('ende.get_value()')
                print(end_node[0].get_value())
            (split_node,List_of_nodes,coverage_of_range)=tree.seach_for_split_node_2(start_node,end_node[0])
            print(coverage_of_range)



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
                crucial_reads.add(node[2])
            else:
                not_crucial_reads[node[2]]=((leaf,node,split_node,List_to_change))
    #crucial is a set of leaf start and end point, correspoinding split node and list of nodes, which have to be changed
    #not _crucial is a dictionary, where for every index the same is stored lik in the crucial set
    return (crucial_reads,not_crucial_reads)

#only for debugging
def detect_crucial_reads_2(BST,max_coverage):
    crucial_reads={}
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
                crucial_reads[node[2]]=((leaf,node,split_node,List_to_change,node[2]))
            else:
                not_crucial_reads[node[2]]=((leaf,node,split_node,List_to_change))
    #crucial is a set of leaf start and end point, correspoinding split node and list of nodes, which have to be changed
    #not _crucial is a dictionary, where for every index the same is stored lik in the crucial set
    return (crucial_reads,not_crucial_reads)



def reduce_readset_via_crucial_and_non_crucial_reads(BST, max_coverage,crucial_set,not_crucial_set):
    '''
    :param: Tree, max_coverage, crucial_set of indices and not_crucial_set, corresponding to reads and  split node
    _:return: First set all former selected crucial indices into the selected readset
    Then go over the leaf list : look if we need to remove something :
    Go over the siblings and look if the index is already in the selected set
    '''
    selected_set=set()
    not_selected_set=set()
    leaf_list=BST.get_leaf_list_of_tree()
    #first add all indices to the selected which are already crucial
    for i in crucial_set:
        selected_set.add(i)
    #assert that all reads are transfered in the selected set
    assert len(selected_set)==len(crucial_set)
    #Go over all reads/intervals

    for leaf in leaf_list:
        new_leaf_coverage=leaf.get_coverage() +leaf.get_balance()
        number_need_to_remove=new_leaf_coverage - max_coverage
        l_siblings=leaf.get_sibling()
        start_nodes= [(sib,i)  for (sib,val,i) in l_siblings if val >leaf.get_value() and i not in selected_set]
        #Store the removed nodes in this list
        just_in_this_case_removed=[]
        while number_need_to_remove>0 and len(start_nodes)!=0:
            (start,index)=start_nodes.pop()
            ((leaf_node,node,split_node,List_to_change))=not_crucial_set[index]
            not_selected_set.add(index)
            just_in_this_case_removed.append(index)
            step_up_balance(List_to_change)
            update_till_root(split_node)
            number_need_to_remove-=1

        if number_need_to_remove>0:
            not_crucial_indices=set(not_crucial_set.keys())
            intersection_between_not_crucial_and_already_selected=not_crucial_indices.intersection(selected_set)
            #Go over this intersecdtion set and look which cover the leaf
            # while number_need_to_remove>0:
            for index in intersection_between_not_crucial_and_already_selected:
                #because end_node is still the whole sibling of the start_node
                # end_node=(BST_leaf_node, value_of_leaf ,index)
                ((start_node,end_node,split_node,List_to_change))=not_crucial_set[index]
                #If leaf lies between start and end node
                if start_node.get_value()<=leaf.get_value() and leaf.get_value()<= end_node[0].get_value():
                    #TODO Not in this case once for the dataset
                    #print('Found interval, where the actual leaf is covered of')
                    not_selected_set.add(index)
                    selected_set.remove(index)
                    step_up_balance(List_to_change)
                    update_till_root(split_node)
                    number_need_to_remove-=1
                if number_need_to_remove==0:
                    break

        if number_need_to_remove==0:
            adding_indices=[i for (sib,val,i) in l_siblings if val > leaf.get_value() and i not in not_selected_set]
            for index in adding_indices:
                selected_set.add(index)
        if number_need_to_remove<0:
            #every read covering this should be added
            adding_indices=[i for (sib,val,i) in l_siblings if val > leaf.get_value() ]
            for index in adding_indices:
                selected_set.add(index)

        new_leaf_coverage=leaf.get_coverage()+leaf.get_balance()
        #Assertion that after reduction the leaf coverage does not exceed the given maximum coverage
        assert new_leaf_coverage<=max_coverage
    return (selected_set,not_selected_set)


#only for debugging
def find_for_this_leaf_all_reads_covering_it(leaf_list,leaf):
    leaf_val=leaf.get_value()
    print('Find for this leaf all reads covering it %d' %leaf_val)

    leaf_inices=leaf.get_index()
    leaf_siblings=leaf.get_sibling()
    for node in leaf_list:
        node_sibling=node.get_sibling()
        node_value=node.get_value()
        if node_value<=leaf_val:
            for (sib,val,i)in node_sibling:
                if val>=leaf_val:
                    print('Found overlapping interval with index %d' %i)
        else:
            #node_value > leaf_val
            for (sib,val,i) in node_sibling:
                if val<= leaf_val:
                    print('Find overlapping interval with index %d' %i)


def step_up_balance(List_to_change):
    '''
    :returns: method which decreases the balance of the given List by one, because the given list includes all reads
    whose balance has to be updatet when this read is pruned
    '''
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
