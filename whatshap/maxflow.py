import math
from whatshap._core import PyRead, PyReadSet
from whatshap.Binary_Search_Tree import Binary_Search_Tree




# TODO
# Implementation of the Interval scheduling problem described in the paper of Veli Mäkinen  "Interval scheduling maximizing minimum coverage "


#IMPORTANT TO NOTICE: NOT SUITABLE FOR PAIRED END READS

#First like in the score_based approach and the random approach: Remove the reads which only cover one variant

#def __init__(self, pruned_readset, max_coverage):
#    self.pruned_readset = PyReadSet()
#    self.max_coverage = max_coverage


#TODO:
#Build up a perfect binary search tree with the delimiters of the intervals(so the reads) as leaves.

#Because we need to keep track of the reads which cover the same regions like the reads discarderd (same problem like
# we had in the priority queue) i implemented the Tree structure by myself with not only parent, siblings, which represent the reads covering this node
#Additionally in the initialization each node stores one value its position and later the attributes are added

def optimize_max_flow_in_BST(readset, BST, max_cov):
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
            print('calling split node')
            #print('In For of optimizing with value of node %d' %leaf_value)
            (split_node, List_to_change) = BST.seach_for_split_node(leaf, end_node)
            print('Called split node with following nodes or reads')
            print('Leaf_node of start %d' % leaf_value)
            print('End node of the read %d' %end_node.get_value())
            if split_node == None:
                print('Found no split node')
            else:
                split_cov = split_node.get_coverage()
                split_balance = split_node.get_balance()
                selection_criterion = BST.is_crucial(split_cov, max_cov, split_balance)
                if selection_criterion:
                    print('Selection Done')
                    pruned_for_ending.append(end_node.get_index())

                #Need to change balance fot the involved nodes
                else:
                    print('Need to remove this read')
                    removed_for_ending.append(end_node.get_index())
                    step_up_balance(BST,List_to_change)
                    #TODO also need to update balance th tree up till the root

#for debugging case return both sets

    return (pruned_for_ending,removed_for_ending)

#after each balance step the siblings should be synchronized again, because of all siblings of a node the balance changed
#TODO : Is that needed? Only problem is because in the split node we get all the different siblings of the same value

def step_up_balance(BST,List_to_change):
    print('Length of the list_to change %d' %len(List_to_change))
    #convert list_to change to set
    new_set=set(List_to_change)
    for l in new_set:
        balance_of_l=l.get_balance()
        l.set_balance(balance_of_l -1)

#Before starting the algorithm we have to build up the binary search tree.
class one_d_range_tree:
    """
	Defines a one dimensional range tree. This has the same structure as a binary search tree, with additional
	properties in the nodes .This is the  minimum and the maximum coverage the balance.

	Nodes ar distinguishedc to leaf nodes and inner nodes

	Problem to deal with are the same coordinates or delimiters in the tree, e.g. Same start variant of a read.
	The idea is to combine the nodes in the way that we have 1 node for all nodes starting at this positions.
	The different reads are then captivated by the sibling attribute which stores the end of the read.
	"""

    def __init__(self, readset):

        '''Build up an list out of the given readset which fullfills the characteristics of a binary search tree '''
        #TODO Call function to construct the tree completely#
        #for each read in readset we need 2 points, a start and an endpoint.

        tree_list = self.build_list(readset)
        full_tree = self.discover_double_and_sibling(tree_list)
        #new_tree=self.develop_layer(full_tree)
        #Need to know if we got uneven number of leaf nodes, after the first layer because then we need to change the structure

        layer_array = []
        complete_tree = self.building_BST_from_leaf_list(0, len(full_tree) - 1, layer_array, full_tree)
        self.complete_tree = complete_tree
        self.leaf_list = full_tree

    def get_complete_tree(self):
        return self.complete_tree

    def get_leaf_list_of_tree(self):
        return self.leaf_list

    #TODO
    def is_crucial(self, read):
        '''
        :param read: Find out if the passed read is crucial
        :return: Value if read is crucial in this  Interval or not an dif it is crucial the reaad is added to the pruned readset
        '''
        return True
        #start_position = read.getposition()
        #end_positions = read.getposition()
        #TODO Look again if it is really ceil
        #if get_min_cov(start_position, end_positions) <= math.ceil(
        #                get_max_cov(start_position, end_positions, self.max_coverage) / 2):
        #    self.pruned_readset.add(read)
        #return False


    def get_min_cov(start, end):
        '''
        returns the minimum coverage in the given interval
        :param start: start position of the interval represents a variant position
        :param end: end position of the interval also equal a variant position
        :return: integer which is the minimum coverage in the whole interval
        '''
        minimum_cov = 0
        return minimum_cov


    def get_max_cov(start, end, max_cov):
        '''
        returns the maximum coverage in the given interval
        :param start: start position of the interval represents a variant position
        :param end: end position of the interval also equal a variant position
        :param max_cov: Given max_cov from whatshap call
        :return: integer which is the maximum coverage in the whole interval, could maximal be the max_cov which was given
         as parameter in the whatshap call
        '''
        maximum_cov = max_cov
        return maximum_cov


    def build_list(self, _ana_readset):
        '''Building up the list of the nodes representing the reads...'''
        list_for_reads = []
        #removing reads which cover less than 2 variants
        indices_of_reads = set(i for i, read in enumerate(_ana_readset) if len(read) >= 2)
        for i in indices_of_reads:
            read_of_index = _ana_readset[i]
            first_pos = read_of_index[0]
            last_pos = read_of_index[len(read_of_index) - 1]
            #initialize the nodes with the position and the position of the sibling
            firs_Node = Leaf_node(first_pos.position, last_pos.position)
            seco_Node = Leaf_node(last_pos.position, first_pos.position)
            list_for_reads.append(firs_Node)
            list_for_reads.append(seco_Node)

        #sort the list for positions or values so that double occuring are in a row and could be removed
        sorted_list = sorted(list_for_reads, key=lambda node: node.value)
        return sorted_list


    def discover_double_and_sibling(self, sorted_list):
        '''
        Detects double occuring leaf nodes erases later and includes the siblings of the erased nodes to the remaining node,
        :param sorted_list: the list of Leaf nodes out of the given readset
        :return: list of leaf nodes with unique leafs and stored siblings and coverage
        '''
        need_to_remove2 = []

        #Not know if coverage count is needed here
        coverage_count = 0

        #Go over the leaf nodes
        iterable = 0
        #print('ITERABLE')
        while iterable != len(sorted_list):
            iter_val = sorted_list[iterable].get_value()
            #coverage_need_to_decrease=0

            #increase coverage if the new node is a start point of some kind of thing
            if (iter_val) < (sorted_list[iterable].get_sibling()[0]):
                coverage_count += 1

            new_var = iterable + 1
            coverage_counter = 0
            while (new_var != len(sorted_list) and iter_val == sorted_list[new_var].get_value() ):

                #If new node start position increase value else remeber to decreas it again.
                if (sorted_list[new_var].get_value()) < (sorted_list[new_var].get_sibling()[0]):
                    coverage_count += 1
                else:
                    coverage_counter += 1
                sorted_list[iterable].add_sibling(sorted_list[new_var].get_sibling())
                if new_var not in need_to_remove2:
                    need_to_remove2.append(new_var)
                new_var += 1

            #Set coverage of the node
            sorted_list[iterable].set_coverage(coverage_count)
            #Reduce counter by the value of the nodes which end in this position
            coverage_count = coverage_count - coverage_counter

            #decrease score if it is the end of the read
            if (iter_val) > (sorted_list[iterable].get_sibling()[0]):
                coverage_count -= 1

            #setting iterable to the point where no doubles occure
            iterable = new_var

        for returned_index in range(0, len(need_to_remove2)):
            to_remove = need_to_remove2.pop()
            sorted_list.pop(to_remove)
        return sorted_list

    def next_layer(self, Going_on, start, tree_list):
        return 0

    def building_BST_from_leaf_list(self, start, end, arr, node_list):
        #print('BST building start %d' %start )
        #print('BST building end %d' %end )

        if (start > end):
            return 0
        if (start == end):
            #root_node=node_list[start]
            return node_list[start]
        else:
            middel = int((start + end) / 2)
            #print('BST Middle %d' %middel)
            (mini, maxi) = self.coverage_of_range(start, end, node_list)
            root_node = BST_node(mini, maxi)
            arr.append(root_node)
            left_node = self.building_BST_from_leaf_list(start, middel, arr, node_list)
            root_node.set_left_child(left_node)
            left_node.set_parent(root_node)
            right_node = self.building_BST_from_leaf_list(middel + 1, end, arr, node_list)
            root_node.set_right_child(right_node)
            right_node.set_parent(root_node)

            return root_node


    def get_all_nodes_though_split_list(self, List_to_change, value_of_sibling):
        length_of_list = len(List_to_change)
        while (length_of_list != 0):
            node = List_to_change.pop()
            node.set_balance(node.get_balance() - 1)
            if node.isLeaf():
                r_child = node.get_right_child()
                l_child = node.get_left_child()
                #when we have a leaf with an bigger value as the end point the balance there has not to be changed.
                if r_child.isLeaf() and r_child.get_value():
                    continue
                l_child.set_balance(node.get_balance() - 1)
                r_child.set_balance(node.get_balance() - 1)
            print('Node :')
            print(node.isLeaf())
            print(node.get_coverage())


    def get_all_nodes_in_the_subtree_without_Those_not_connected_to_sibling(self, split_node_of_read, value_of_sibling,
                                                                            List_to_change):
        print('One_call_of_the_algorithm')
        balance = split_node_of_read.get_balance()
        split_node_of_read.set_balance(balance - 1)
        Right_child_Leaf = False
        Left_child_Leaf = False
        # while not(Right_child_Leaf and Left_child_Leaf):
        #    l_child=split_node_of_read.get_left_child()
        #   r_child=split_node_of_read.get_right_child()
        #  if l_child.isLeaf():
        #     Left_child_Leaf=True
        #     continue
        # if r_child.isLeaf():
        #     Right_child_Leaf=True
        #     continue
        #Right_child_Leaf=True
        #Left_child_Leaf=True


        List_of_nodes_to_change = []
        for i in List_to_change:
            print('For each node coverage, parent, rchild, and left child')
            print(i.get_coverage())
            print(i.get_parent().get_coverage())
            print(i.get_right_child().get_coverage())
            print(i.get_left_child().get_coverage())
            if i.get_right_child().isLeaf():
                print(i.get_right_child().get_value())
            if i.get_left_child().isLeaf():
                print(i.get_left_child().get_value())
        return None


    #TODO does nothing in the moment
    def balance_control(self, split_node_of_read, List_to_change, l_node, value_of_sibling):
        '''
        Should update the balance in the tree, because the read starting at l_node and ending on value of sibling is pruned
        out of the original readset
        :param split_node_of_read: Found split node, so we have to check from this point on up an d down
        :param List_to_change: List which stores which inner BST nodes are followed to get to the split node
        :param l_node: start node of the removed read
        :param value_of_sibling:  value of the end node of the read, BUT NOT THE NODE ITSELF
        :return: Nothing a changed tree with changed balance
        '''
        #List_of_nodes=self.get_all_nodes_in_the_subtree_without_Those_not_connected_to_sibling(split_node_of_read,value_of_sibling,List_to_change)
        self.get_all_nodes_though_split_list(List_to_change, value_of_sibling)
        return None


    #TODO it should be efficienter to store immediatly all nodes by the search for the split node

    def get_split_node(self, start_node, value_of_sibling):
        '''returns the split node, and a  list of nodes which coverage or balance maybe has to be updated, when the read is deleted'''
        List_of_Leafs = []
        List_of_nodes_to_change_by_removal = []
        Found_split_node = False
        split_node = None
        while not Found_split_node:
            List_of_Leafs = []
            p_node = start_node.get_parent()
            List_of_nodes_to_change_by_removal.append(p_node)
            r_child = p_node.get_right_child()
            Found_leaf_list = r_child.get_Leaf_nodes_of_subtree(List_of_Leafs)
            Found_leaf_list_values = [leaf.get_value() for leaf in Found_leaf_list]
            if value_of_sibling in Found_leaf_list_values:
                #print('In if')
                Found_split_node = True
                split_node = p_node
            else:
                #print('In else')
                start_node = p_node
                #print('Hi')
                #print('Start Node coverage: ')
                #print(start_node.get_coverage())
                #To get only one iteration
                #Found_split_node=True
        return (split_node, List_of_nodes_to_change_by_removal)

    def reducing_readset_for_max_coverage(self, reads, max_cov):
        '''

        :param reads: Original readset
        :param max_cov: maximum coverage of the
        :return:
        '''
        pruned_readset = reads
        Need_to_remove = []
        #tree=self.get_complete_tree()
        root_node = self.get_complete_tree()
        leaf_list_of_tree = self.get_leaf_list_of_tree()
        #go over leafs =reads
        for i, l_node in enumerate(leaf_list_of_tree):
            leaf_value = l_node.get_value()
            siblings = l_node.get_sibling()
            #print('L_node . so read node coverage %d' %l_node.get_coverage())
            #print('L_node . value %d' %leaf_value)
            #   new_sibling=sorted(siblings)
            #look only at the start points
            only_end_points = [sib for sib in sorted(siblings) if sib > leaf_value]
            for end_part_of_read in only_end_points:
                #print('Leaf node%d'% end_part_of_read)
                value_of_sibling = end_part_of_read
                print('Calling split node ')
                (split_node_of_read, List_to_change) = self.get_split_node(l_node, value_of_sibling)
                #print('Split node')
                #print(split_node_of_read)
                if split_node_of_read != None:
                    #min and max coverage as combination of coverage with the balance
                    min_coverage_split = split_node_of_read.get_min_coverage() + split_node_of_read.get_balance()
                    max_coverage_split = split_node_of_read.get_max_coverage() + split_node_of_read.get_balance()
                    if min_coverage_split > max_cov:
                        Need_to_remove.append((l_node, value_of_sibling))
                        #TODO Here to add the balance of all nodes which are involved
                        self.balance_control(split_node_of_read, List_to_change, l_node, value_of_sibling)
                        print('In If')
                    else:
                        print('In ELSE')
                        if max_coverage_split > max_cov:
                            read = (l_node, value_of_sibling)
                            #TODO not working also
                            self.is_crucial(read)

                            #print('Split node coverage')
                            #print(split_node_of_read.get_coverage())
                            #print('len of List to change by removal%d' %len(List_to_change))
                            #for i in List_to_change:
                            #print('List to change')
                            #print(i.get_coverage())
                            #print(i.get_balance())
                            #print(only_end_points)
                            #print('Only start points')
                            #print(l_node.get_parent().get_coverage())
        return pruned_readset

    def optimize_max_flow_in_BST(readset, BST, max_cov):
        pruned_readset = []
        leaf_list = BST.get_leaf_list_of_tree()
        for i, leaf in enumerate(leaf_list):
            leaf_value = leaf.get_value()
            siblings = leaf.get_sibling()
            only_end_points = [sib for sib in sorted(siblings) if sib.get_value() > leaf_value]
            BST.synchronize_sibling_with_same_value(only_end_points)


    #Probably not needed
    def develop_layer(self, node_list):
        '''

        :param node_list:
        :return:
        '''
        last_layer = len(node_list)
        print('last layer  %d' % last_layer)
        parent_node_list = []
        parent_node_list = []

        #for i in range(0,len(node_list),2):
        #    print('In Range')
        #    print(i)
        #    print(node_list[i].get_value())
        #    print(node_list[i+1].get_coverage())

        i = 0
        runvariable = 0

        while i < last_layer:
            print('WHILE %d' % i)
            first_node = node_list[i]
            second_node = node_list[i + 1]
            #coverage_1=first_node.get_coverage()
            parent_node = Inner_node(first_node, second_node, first_node.get_coverage(), second_node.get_coverage())
            first_node.set_leaf_node_attributes(last_layer + runvariable)
            second_node.set_leaf_node_attributes(last_layer + runvariable)

            parent_node_list.append(parent_node)
            i += 2
            #runvariable represents the index in the second list  which is later concatenated
            runvariable += 1

        for parent in parent_node_list:
            print(parent.get_min_cov())
            print(parent.get_max_cov())
        print('Leafnodes')
        for leafnodes in node_list:
            print(leafnodes.get_value())
            print(leafnodes.get_parent())

        print('NEW LAYER LIST ')
        new_layer_list = node_list + parent_node_list
        p = 0
        for k in new_layer_list:
            print('K is leaf %d' % k.isLeaf())
            if k.isLeaf():
                print(p)
                print(k.get_value())
            else:
                print(p)
                print(k.get_min_cov())
            p += 1
        print(new_layer_list)

        return node_list


    def coverage_of_range(self, start, stop, nodelist):
        #initializing max and min coverage by -1 because the coverage could not be negativebut minimum has to be high
        maximum = -1
        minimum = 50
        #go over the nodelist
        #Need to add 1 because in the loop the last element is not considered
        for i in range(start, stop + 1):

            #if Leaf nodes
            if (nodelist[i].isLeaf()):
                coverage = nodelist[i].get_coverage()
                if coverage > maximum:
                    maximum = coverage
                if coverage < minimum:
                    minimum = coverage
            #if inner nodes, we have 2 coverages
            else:
                (min_cov, max_cov) = nodelist[i].get_coverage()
                if min_cov < minimum:
                    minimum = min_cov
                if max_cov > maximum:
                    maximum = max_cov
        return (minimum, maximum)


#Nodes represent Nodes in the binary search tree.
#Differentiate between inner nodes and leaf nodes


class BST_node:
    def __init__(self, minimum, maximum):
        self.balance = 0
        self.min_coverage = minimum
        self.max_coverage = maximum

    def get_coverage(self):
        return (self.min_coverage, self.max_coverage)

    def set_left_child(self, left_Node):
        self.left_child = left_Node

    def set_right_child(self, right_Node):
        self.right_child = right_Node

    def get_left_child(self):
        return self.left_child

    def get_right_child(self):
        return self.right_child

    def get_balance(self):
        return self.balance


    def set_balance(self, new_balance):
        self.balance = new_balance

    def get_min_coverage(self):
        return self.min_coverage

    def get_max_coverage(self):
        return self.max_coverage

    def set_parent(self, node):
        self.parent = node

    def get_parent(self):
        return self.parent

    def get_all_leaf_nodes_of_subtree(self):
        while not self.isLeaf():
            print('In While loop')
        #            node=self.get_right_child()
        #return self.get_coverage()
        #Leaf_list=[]
        #while not self.isLeaf():
        #    rc=self.get_right_child()
        #    lc=self.get_left_child()
        #    list_rc=rc.get_all_leaf_nodes_of_subtree()
        #    list_lc=lc.get_all_leaf_nodes_of_subtree()
        #    Leaf_list.append(list_lc)
        #    Leaf_list.append(list_rc)
        #if self.isLeaf():
        #    Leaf_list.append(self)
        #return Leaf_list
        return 0

    def isLeaf(self):
        return False

    def return_node(self):

        return self

    def get_Leaf_nodes_of_subtree(self, List_of_Leafs):
        '''
        :param List_of_Leafs: Gets an List of Leafs which are already discovered
        :return:List of Leafes of the childs of this node
        '''
        r_child = self.get_right_child()
        l_child = self.get_left_child()
        if r_child.isLeaf():
            List_of_Leafs.append(r_child)
        else:
            r_child.get_Leaf_nodes_of_subtree(List_of_Leafs)
        if l_child.isLeaf():
            List_of_Leafs.append(l_child)
        else:
            l_child.get_Leaf_nodes_of_subtree(List_of_Leafs)
        return List_of_Leafs


class Leaf_node:
    def __init__(self, value, sibling):
        '''
        :param value: Corresponds to SNP either first or last position in a read
        :param sibling: First only the node which represents the other delimiter of the included read,
                if Leaf nodes with same positions exits they are combined and then the sibling changes to a list
        :return: A leaf node with at attributes value and at least a sibling.
        '''
        self.value = value
        adding_sibling = []
        adding_sibling.append(sibling)
        self.sibling = adding_sibling

    def set_parent(self, parent):
        '''Setting the node attributes
        Coverage is exactly the number of sibling the node has
        :param parent: The index of the parent node in the array
        :return
        '''
        self.balance = 0
        self.parent = parent

    #Does not work with len(self.sibling) because then the coverage from other reads which do not start or end there is neglected,

    def set_coverage(self, coverage):
        self.coverage = coverage


    def get_parent(self):
        return self.parent

    def get_coverage(self):
        return self.coverage

    def get_balance(self):
        return self.balance

    def get_value(self):
        return self.value

    def get_sibling(self):
        #returns list of sibling, either one ore more if the position of the node occured more than once
        return self.sibling

    def add_sibling(self, value_list):
        adding_sibling = self.sibling
        for i in value_list:
            adding_sibling.append(i)
        self.sibling = adding_sibling

    def isLeaf(self):
        return True

    def get_all_leaf_nodes_of_subtree(self):
        empty_list = []
        return empty_list


    def get_Leaf_nodes_of_subtree(self, List_of_Leafs):
        '''
        :param List_of_Leafs: Gets an List of Leafs which are already discovered
        :return:List of Leafes of the childs of this node
        '''
        List_of_Leafs.append(self)
        return List_of_Leafs

