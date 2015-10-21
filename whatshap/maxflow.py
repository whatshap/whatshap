import math
from whatshap._core import PyRead,PyReadSet

# TODO
#Implementation of the Interval scheduling problem described in the paper of Veli Mäkinen  "Interval scheduling maximizing minimum coverage "


#IMPORTANT TO NOTICE: NOT SUITABLE FOR PAIRED END READS

#First like in the score_based approach and the random approach: Remove the reads which only cover one variant

def __init__(self, pruned_readset, max_coverage):
    self.pruned_readset = PyReadSet()
    self.max_coverage = max_coverage


def is_crucial(self, read):
    '''
    :param read: Find out if the passed read is crucial
    :return: Value if read is crucial in this  Interval or not an dif it is crucial the reaad is added to the pruned readset
    '''
    start_position = read.getposition()
    end_positions = read.getposition()
    #TODO Look again if it is really ceil
    if min_cov(start_position, end_positions) <= math.ceil(
                    max_cov(start_position, end_positions, self.max_coverage) / 2):
        self.pruned_readset.add(read)
    return False


def min_cov(start, end):
    '''
    returns the minimum coverage in the given interval
    :param start: start position of the interval represents a variant position
    :param end: end position of the interval also equal a variant position
    :return: integer which is the minimum coverage in the whole interval
    '''
    minimum_cov = 0
    return minimum_cov


def max_cov(start, end, max_cov):
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
        print('tree_list')
        for i in range(0,len(full_tree)):
            val=tree_list[i].get_value()
            sib=tree_list[i].get_sibling()
            print('val & sib')
            print(val)
            print(sib)



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
            print('Need to look if the positions in the sibling object of the read change with the removing of the double occurences should be no problem because we store positions and not indices ')
            print('First Node : own position  and sibling position')
            print(first_pos.position)
            print(last_pos.position)
            list_for_reads.append(firs_Node)
            list_for_reads.append(seco_Node)

        #sort the list for positions or values so that double occuring are in a row and could be removed
        sorted_list=sorted(list_for_reads,key=lambda node :node.value)
        return sorted_list


    def discover_double_and_sibling(self,sorted_list):
        '''
        Detects double occuring leaf nodes erases later and includes the siblings of the erased nodes to the remaining node,
        :param sorted_list: the list of Leaf nodes out of the given readset
        :return: list of leaf nodes with unique leafs and stored siblings
        '''
        need_to_remove=[]
        need_to_remove2=[]

        #Not know if coverage count is needed here
        coverage_count=0
        #TODO : Not known if needed
        min_cov_in_the_range= 0
        max_cov_in_the_range=0
        unique_node_val= (0,0)
        Going_on = True

        #Go over the leaf nodes

        ##TRY ALTERNATIVe :
        iterable=0
        print('ITERABLE')
        while iterable !=len(sorted_list):
            print('In WHILE Loop')
            print(sorted_list[iterable].get_value())
            iter_val= sorted_list[iterable].get_value()
            coverage_need_to_decrease=0

            #increase coverage if the new node is a start point of some kind of thing
            if (iter_val)<(sorted_list[iterable].get_sibling()[0]):
                coverage_count +=1
            else:
                coverage_need_to_decrease+=1

            new_var= iterable+1

            while (new_var!= len(sorted_list) and iter_val==sorted_list[new_var].get_value() ):
                print('In WHIle LOOP WHILE DOUNF DOUBLE ')
                new_var +=1
                #start
                #If new node start position increase value else remeber to decreas it again.
                if (sorted_list[new_var].get_value())<(sorted_list[new_var].get_sibling()[0]):
                    coverage_count +=1
                else:
                    coverage_need_to_decrease +=1

                sorted_list[iterable].add_sibling(sorted_list[new_var].get_sibling())
                if new_var not in need_to_remove2:
                    need_to_remove2.append(new_var)
            print('Found folloing coverage')
            print(coverage_count)
            coverage_count=coverage_count-coverage_need_to_decrease


            #setting iterable to the point where no doubles occure
            iterable=new_var
        print()

        for returned_index in range(0,len(need_to_remove2)):
            to_remove=need_to_remove2.pop()
            sorted_list.pop(to_remove)
        #print('Print i s siblings' )
        #for i in range(0,len(sorted_list)):
        #    print(sorted_list[i].get_sibling())
        return sorted_list


Without_Bridging


'''
        for i in range(0,len(sorted_list)):
            print('len(sorted_list)')
            print(len(sorted_list))
            #print(sorted_list[len(sorted_list)])
            print('I in first for loop')
            print(i)
            print(coverage_count)
            #first position Node :
            actual_val= sorted_list[i].get_value()


            if actual_val !=(unique_node_val[0]):
                print('In If unique node val')
                print(unique_node_val)

            #so the actual_val represents a start position
            if (actual_val)<(sorted_list[i].get_sibling()[0]):
                coverage_count +=1
            print('After adding coverage')
            print(coverage_count)
            unique_node_val= (actual_val,coverage_count)
            #define a new Node
            j=i+1
            #while nodes with same postion exists, then we need to remove them from the list of leaves and store the sibling
            #in a list
            if j !=len(sorted_list):
                while (sorted_list[j].get_value()==actual_val):
                    unique_node_val=(actual_val,coverage_count+1)
                    #improve COVerage
                    #TODO Try here like on the block !!!!

                    coverage_count+=1
                    #print('Append j at the moment')
                    #print(j)
                    #need_to_remove.append(j)
                    sorted_list[i].add_sibling(sorted_list[j].get_sibling())
                    if j not in need_to_remove:
                        #print('In if case')
                        need_to_remove.append(j)
                    j+=1
                    #print('I and J in the while loop')
                    #print(i)
                    #print(j)
                    #print('New j ')
                    #print(j)#


            #means you are at the right border of the read so the coverage has to be decrease
            if (actual_val)>(sorted_list[i].get_sibling()[0] and sorted_list[i+1].get_value()!= actual_val):
                coverage_count -=1
            print('After removing cov')
            print(coverage_count)


            #if no double occurences then add 1 otherwise skip the one double occuring
            #i=j
            #print('new i')
            #print(i)
        print('Last unique node')
        print(unique_node_val)
        #print('NEED TO REMOVE IT ')
        #print(need_to_remove)

        for returned_index in range(0,len(need_to_remove)):
            #print('returned_index')
            to_remove=need_to_remove.pop()
            #sorted_list[to_remove]
            #need_to_remove.remove(sorted_list[to_remove])
            #print(need_to_remove)

            #print(to_remove)
            #print(need_to_remove)
            sorted_list.pop(to_remove)

        #print('Print i s siblings' )
        #for i in range(0,len(sorted_list)):
        #    print(sorted_list[i].get_sibling())

'''





#Nodes represent Nodes in the binary search tree.
#Differentiate between inner nodes and leaf nodes

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

    def set_leaf_node_attributes(self, parent,coverage):
        '''Setting the node attributes
        Coverage is exactly the number of sibling the node has
        :param parent: The index of the parent node in the array
        :return
        '''
        self.coverage =coverage
        #Does not work with len(self.sibling) because then the coverage from other reads which do not start or end there is neglected,
        self.balance=0
        self.parent=parent


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
        return self.max_cov()

    def get_balance(self):
        return self.balance

    def parent(self):
        return self.parent

    def get_left_child(self):
        return self.left

    def get_right_child(self):
        self.right

