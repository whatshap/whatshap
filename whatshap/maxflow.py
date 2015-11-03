import math
from whatshap._core import PyRead,PyReadSet

# TODO
#Implementation of the Interval scheduling problem described in the paper of Veli Mäkinen  "Interval scheduling maximizing minimum coverage "


#IMPORTANT TO NOTICE: NOT SUITABLE FOR PAIRED END READS

#First like in the score_based approach and the random approach: Remove the reads which only cover one variant

#def __init__(self, pruned_readset, max_coverage):
#    self.pruned_readset = PyReadSet()
#    self.max_coverage = max_coverage


def is_crucial(self, read):
    '''
    :param read: Find out if the passed read is crucial
    :return: Value if read is crucial in this  Interval or not an dif it is crucial the reaad is added to the pruned readset
    '''
    start_position = read.getposition()
    end_positions = read.getposition()
    #TODO Look again if it is really ceil
    if get_min_cov(start_position, end_positions) <= math.ceil(
                    get_max_cov(start_position, end_positions, self.max_coverage) / 2):
        self.pruned_readset.add(read)
    return False


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


#TODO:
#Build up a perfect binary search tree with the delimiters of the intervals(so the reads) as leaves.

#Because we need to keep track of the reads which cover the same regions like the reads discarderd (same problem like
# we had in the priority queue) i implemented the Tree structure by myself with not only parent, siblings, which represent the reads covering this node
#Additionally in the initialization each node stores one value its position and later the attributes are added



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

        tree_list=self.build_list(readset)
        full_tree=self.discover_double_and_sibling(tree_list)
        #new_tree=self.develop_layer(full_tree)
        #Need to know if we got uneven number of leaf nodes, after the first layer because then we need to change the structure

        layer_array=[]
        complete_tree=self.building_BST_from_leaf_list(0,len(full_tree)-1,layer_array,full_tree)
        self.complete_tree=complete_tree
        self.leaf_list=full_tree

    def get_complete_tree(self):
        return self.complete_tree

    def get_leaf_list_of_tree(self):
        return self.leaf_list

    def build_list(self,_ana_readset):
        '''Building up the list of the nodes representing the reads...'''
        list_for_reads=[]
        #removing reads which cover less than 2 variants
        indices_of_reads = set(i for i, read in enumerate(_ana_readset) if len(read) >= 2)
        for i in indices_of_reads:

            read_of_index= _ana_readset[i]
            first_pos=read_of_index[0]
            last_pos=read_of_index[len(read_of_index)-1]
            #initialize the nodes with the position and the position of the sibling
            firs_Node=Leaf_node(first_pos.position,last_pos.position)
            seco_Node= Leaf_node(last_pos.position,first_pos.position)
            list_for_reads.append(firs_Node)
            list_for_reads.append(seco_Node)

        #sort the list for positions or values so that double occuring are in a row and could be removed
        sorted_list=sorted(list_for_reads,key=lambda node :node.value)
        return sorted_list


    def discover_double_and_sibling(self,sorted_list):
        '''
        Detects double occuring leaf nodes erases later and includes the siblings of the erased nodes to the remaining node,
        :param sorted_list: the list of Leaf nodes out of the given readset
        :return: list of leaf nodes with unique leafs and stored siblings and coverage
        '''
        need_to_remove2=[]

        #Not know if coverage count is needed here
        coverage_count=0

        #Go over the leaf nodes
        iterable=0
        #print('ITERABLE')
        while iterable !=len(sorted_list):
            iter_val= sorted_list[iterable].get_value()
            #coverage_need_to_decrease=0

            #increase coverage if the new node is a start point of some kind of thing
            if (iter_val)<(sorted_list[iterable].get_sibling()[0]):
                coverage_count +=1

            new_var= iterable+1
            coverage_counter= 0
            while (new_var!= len(sorted_list) and iter_val==sorted_list[new_var].get_value() ):

                #If new node start position increase value else remeber to decreas it again.
                if (sorted_list[new_var].get_value())<(sorted_list[new_var].get_sibling()[0]):
                    coverage_count +=1
                else:
                    coverage_counter+=1
                sorted_list[iterable].add_sibling(sorted_list[new_var].get_sibling())
                if new_var not in need_to_remove2:
                    need_to_remove2.append(new_var)
                new_var +=1

            #Set coverage of the node
            sorted_list[iterable].set_coverage(coverage_count)
            #Reduce counter by the value of the nodes which end in this position
            coverage_count=coverage_count-coverage_counter

            #decrease score if it is the end of the read
            if (iter_val)>(sorted_list[iterable].get_sibling()[0]):
                coverage_count -=1

            #setting iterable to the point where no doubles occure
            iterable=new_var

        for returned_index in range(0,len(need_to_remove2)):
            to_remove=need_to_remove2.pop()
            sorted_list.pop(to_remove)
        return sorted_list

    def next_layer(self,Going_on,start,tree_list):
        return 0

    def building_BST_from_leaf_list(self,start,end,arr,node_list):
        #print('BST building start %d' %start )
        #print('BST building end %d' %end )

        if (start>end):
            return 0
        if (start==end):
            #root_node=node_list[start]
            return node_list[start]
        else:
            middel= int((start+end)/2)
            #print('BST Middle %d' %middel)
            (mini,maxi)=self.coverage_of_range(start,end,node_list)
            root_node=BST_node(mini,maxi)
            arr.append(root_node)
            left_node=self.building_BST_from_leaf_list(start,middel,arr,node_list)
            root_node.set_left_child(left_node)
            left_node.set_parent(root_node)
            right_node=self.building_BST_from_leaf_list(middel+1,end,arr,node_list)
            root_node.set_right_child(right_node)
            right_node.set_parent(root_node)

            return root_node



    def get_split_node(self,start_node,value_of_sibling):
        '''returns the split node, and a  list of nodes which coverage or balance maybe has to be updated, when the read is deleted'''
        List_of_Leafs=[]
        List_of_nodes_to_change_by_removal=[]
        Found_split_node=False
        split_node= None
        while not Found_split_node:
            List_of_Leafs=[]
            p_node=start_node.get_parent()
            List_of_nodes_to_change_by_removal.append(p_node)
            r_child=p_node.get_right_child()
            Found_leaf_list=r_child.get_Leaf_nodes_of_subtree(List_of_Leafs)
            Found_leaf_list_values=[leaf.get_value() for leaf in Found_leaf_list]
            if value_of_sibling in Found_leaf_list_values:
                print('In if')
                Found_split_node=True
                split_node=p_node
            else:
                print('In else')
                start_node=p_node
            print('Hi')
            print('Start Node coverage: ')
            print(start_node.get_coverage())
            #To get only one iteration
            #Found_split_node=True
        return (split_node,List_of_nodes_to_change_by_removal)






    #Probably not needed
    def develop_layer(self,node_list):
        '''

        :param node_list:
        :return:
        '''
        last_layer=len(node_list)
        print('last layer  %d' %last_layer)
        parent_node_list=[]
        parent_node_list=[]

        #for i in range(0,len(node_list),2):
        #    print('In Range')
        #    print(i)
        #    print(node_list[i].get_value())
        #    print(node_list[i+1].get_coverage())

        i=0
        runvariable=0

        while i<last_layer:
            print('WHILE %d'% i )
            first_node= node_list[i]
            second_node=node_list[i+1]
            #coverage_1=first_node.get_coverage()
            parent_node=Inner_node(first_node,second_node,first_node.get_coverage(),second_node.get_coverage())
            first_node.set_leaf_node_attributes(last_layer+runvariable)
            second_node.set_leaf_node_attributes(last_layer+runvariable)

            parent_node_list.append(parent_node)
            i+=2
            #runvariable represents the index in the second list  which is later concatenated
            runvariable+=1

        for parent in parent_node_list:
            print(parent.get_min_cov())
            print(parent.get_max_cov())
        print('Leafnodes')
        for leafnodes in node_list:
            print(leafnodes.get_value())
            print(leafnodes.get_parent())

        print('NEW LAYER LIST ')
        new_layer_list=node_list + parent_node_list
        p=0
        for k in new_layer_list:
            print('K is leaf %d' %k.isLeaf())
            if k.isLeaf():
                print(p)
                print(k.get_value())
            else:
                print(p)
                print(k.get_min_cov())
            p +=1
        print(new_layer_list)

        return node_list


    def coverage_of_range(self,start,stop,nodelist):
        #initializing max and min coverage by -1 because the coverage could not be negativebut minimum has to be high
        maximum = -1
        minimum = 50
        #go over the nodelist
        #Need to add 1 because in the loop the last element is not considered
        for i in range(start,stop+1):

            #if Leaf nodes
            if (nodelist[i].isLeaf()):
                coverage=nodelist[i].get_coverage()
                if coverage>maximum:
                    maximum=coverage
                if coverage<minimum:
                    minimum=coverage
            #if inner nodes, we have 2 coverages
            else:
                (min_cov,max_cov)=nodelist[i].get_coverage()
                if min_cov<minimum :
                    minimum=min_cov
                if max_cov>maximum:
                    maximum=max_cov
        return (minimum,maximum)

#Nodes represent Nodes in the binary search tree.
#Differentiate between inner nodes and leaf nodes


class BST_node:
    def __init__(self,minimum,maximum):
        self.balance=0
        self.min_coverage=minimum
        self.max_coverage=maximum

    def get_coverage(self):
        return (self.min_coverage,self.max_coverage)

    def set_left_child(self,left_Node):
        self.left_child=left_Node

    def set_right_child(self,right_Node):
        self.right_child=right_Node

    def get_left_child(self):
        return self.left_child

    def get_right_child(self):
        return self.right_child

    def get_balance(self):
        return self.balance

    def get_min_coverage(self):
        return self.min_coverage

    def get_max_coverage(self):
        return self.max_coverage

    def set_parent(self,node):
        self.parent= node

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

    def get_Leaf_nodes_of_subtree(self,List_of_Leafs):
        '''
        :param List_of_Leafs: Gets an List of Leafs which are already discovered
        :return:List of Leafes of the childs of this node
        '''
        r_child=self.get_right_child()
        l_child=self.get_left_child()
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
    def __init__(self,value, sibling):
        '''
        :param value: Corresponds to SNP either first or last position in a read
        :param sibling: First only the node which represents the other delimiter of the included read,
                if Leaf nodes with same positions exits they are combined and then the sibling changes to a list
        :return: A leaf node with at attributes value and at least a sibling.
        '''
        self.value=value
        adding_sibling=[]
        adding_sibling.append(sibling)
        self.sibling=adding_sibling

    def set_parent(self, parent):
        '''Setting the node attributes
        Coverage is exactly the number of sibling the node has
        :param parent: The index of the parent node in the array
        :return
        '''
        self.balance=0
        self.parent=parent

#Does not work with len(self.sibling) because then the coverage from other reads which do not start or end there is neglected,

    def set_coverage(self,coverage):
        self.coverage=coverage


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

    def add_sibling(self,value_list):
        adding_sibling= self.sibling
        for i in value_list:
            adding_sibling.append(i)
        self.sibling = adding_sibling

    def isLeaf(self):
        return True

    def get_all_leaf_nodes_of_subtree(self):
        empty_list=[]
        return empty_list



    def get_Leaf_nodes_of_subtree(self,List_of_Leafs):
        '''
        :param List_of_Leafs: Gets an List of Leafs which are already discovered
        :return:List of Leafes of the childs of this node
        '''
        List_of_Leafs.append(self)
        return List_of_Leafs

class Inner_node:
    def __init__(self,left_child,right_child,cov_1,cov_2):
        '''initializing a inner node.
        difference from the leaf is the missing position value and the included child information
        the parent information will be defined later, if it is needed.
        In the initialization process the balance is set to 0

        :param left_child: Index of the left child in the tree
        :param right_child: Indey of the right child in the tree
        :param cov_1: Coverage of one child
        :param cov_2: Coverage of the second child
        :return: a inner node with different attributes
        '''
        self.left=left_child
        self.right=right_child
        self.min_cov=min(cov_1,cov_2)
        self.max_cov=max(cov_1,cov_2)
        self.balance=0


    def set_parent(self, parent):
        self.parent=parent


    def get_min_cov(self):
        return self.min_cov

    def get_max_cov(self):
        return self.max_cov

    def get_balance(self):
        return self.balance

    def get_parent(self):
        return self.parent

    def get_left_child(self):
        return self.left

    def get_right_child(self):
        self.right

    def get_coverage(self):
        return (min_cov , max_cov)

    def isLeaf(self):
        return False

