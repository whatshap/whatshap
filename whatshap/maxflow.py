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

#Because we could not start from first position to last position or delimiter in the genome like suggested
#start with the first read, and go over the reads, because the reads are already sorted by their position in the genome



#Because we need to keep track of the reads which cover the same regions like the reads discarderd (same problem like
# we had in the priority queue) i implemented the Tree structure by myself with not only parent, left and right child
# but also left and right_sibling, which could be detected by a indix in a list or array structure,,,,
#Additionally in the initialization each node stores only one value and later it stores 2 values..

#Also additionally a balance counter...



class one_d_range_tree:
    """
	Defines a one dimensional range tree. This has the same structure as a binary search tree, with additional
	properties in the nodes .This is the  minimum and the maximum coverage the balance.

    Need additional the Coverage Monitor to get the true min and max values of an interval, because the crucial poitn could
    lie between the delimiter, or the variant positions.

	Each node has a parent, left and right child and  a sibling direct.

	Problem to deal with are the same coordinates or delimiters in the tree, e.g. Same start point of a read.
	I will just ignore them because the coice which reads are removed is random and we just have to filter out the crucial ones first.

	"""

    def __init__(self, readset):

        '''Build up an array out of the given readset'''
    #Call function to construct the tree completely#
        #for each read in readset we need 2 points, a start and an endpoint.
        coverage_monitor= Min_Max_CovMonitor(readset.get_positions())

        #self.nodes = {x: Node(x, None) for x in readset}
        first_list= self.build_list(readset)
        print('FIRST  LIST ')
        new_tree= self.discover_double_add_coverage_and_parent_or_sibling(first_list)


    def build_list(_ana_readset):
        '''Building up the list of the nodes representing the reads...'''
        list_for_reads=[]
        #readset.get_positions()
        indices_of_reads = set(i for i, read in enumerate(_ana_readset) if len(read) >= 2)
        for i in indices_of_reads:
            read_of_index= _ana_readset[i]
            #print('read_of_index')
            #print(read_of_index)
            first_pos=read_of_index[0]
            last_pos=read_of_index[len(read_of_index)-1]
            #print(len(read_of_index))
            #print(first_pos.position)
            #print(last_pos.position)
            firs_Node=Node(first_pos.position)
            seco_Node= Node(last_pos.position)
            list_for_reads.append(firs_Node)
            list_for_reads.append(seco_Node)

        #print('List of Reads in Structure')
        #print(list_for_reads)


        #for nodes in list_for_reads:
        #    print(nodes.get_value())

        sorted_list=sorted(list_for_reads,key=lambda node :node.value)

        #print('After sorting')

        #for nodes in sorted_list:
        #    print(nodes.get_value())

        return sorted_list

        #Done in the init of the nodes
        #discover_double_add_coverage_and_parent_or_sibling(sorted_list)

        #minimum_val=0
        #coverage=0
        #for read in readset:
            #coverage +=1
        #    print(read.getVariantCount())
            #list_for_reads.append(Node(read.firstPosition()))
            #print('First Positions')
            #print(read.firstPosition())
            #print('Last Position')
            #print(read.lastPosition())
        #    list_for_reads.append(Node(read.lastPosition()))
            #coverage +=1

        #list_for_reads.sort()

    def discover_double_add_coverage_and_parent_or_sibling(self,sorted_list):
        for i in range(0,len(sorted_list)-1):
            actual_val= sorted_list[i].get_value()
            j=i+1
            while (sorted_list[j]==actual_val):
                print('Found double ')
            i +=1







class Node:
    def __init__(self,value):
        self.value=value

 #   def __init__(self, min_cov, max_cover, parent, sibling, balance):
 #       self.min_cov = min_cov
 #       self.sibling=sibling
 #       self.balance=balance
 #       self.max_cov=max_cover
 #       self.parent=parent

    def Node_attributes(self,coverage, parent, sibling, balance):
        self.min_cov = min_cov
        self.sibling=sibling
        self.balance=balance
        self.parent=parent


    def get_min_cov(self):
        return self.min_cov

    #def get_max_cov(self):
     #   return self.max_cov()

    def get_balance(self):
        return self.balance

    def get_value(self):
        return self.value


#class Min_Max_CovMonitor:
#    def __init__(self, length):
#        self.coverage = [0] * length
#
#    def max_coverage_in_range(self, begin, end):
#        return max(self.coverage[begin:end])
#
#    def min_coverage_in_range(self,begin,end):
#        return min(self.coverage[begin:end])
#
#    def add_read(self, begin, end):
#        for i in range(begin, end):
#            self.coverage[i] += 1
#
#    def remove_read(self,begin,end):
#        for i in range(begin,end):
#            self.coverage[i] -=1


