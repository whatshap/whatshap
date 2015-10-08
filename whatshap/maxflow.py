import math
from ..core import Read, ReadSet

#TODO
#Implementation of the Interval scheduling problem described in the paper of Veli Mäkinen  "Interval scheduling maximizing minimum coverage "


#First like in the score_based approach and the random approach: Remove the reads which only cover one variant

def __init__(self, pruned_readset,max_coverage):
    self.pruned_readset=ReadSet()
    self.max_coverage= max_coverage



def is_crucial(self,read):
    '''
    :param read: Find out if the passed read is crucial
    :return: Value if read is crucial in this  Interval or not an dif it is crucial the reaad is added to the pruned readset
    '''
    start_position= read.getposition()
    end_positions=read.getposition()
    #TODO Look again if it is really ceil
    if min_cov(start_position,end_positions)<= math.ceil(max_cov(start_position,end_positions,self.max_coverage) /2):
        self.pruned_readset.add(read)
    return False




def min_cov(start,end):
    '''
    returns the minimum coverage in the given interval
    :param start: start position of the interval represents a variant position
    :param end: end position of the interval also equal a variant position
    :return: integer which is the minimum coverage in the whole interval
    '''
    minimum_cov=0
    return minimum_cov


def max_cov(start,end,max_cov):
    '''
    returns the maximum coverage in the given interval
    :param start: start position of the interval represents a variant position
    :param end: end position of the interval also equal a variant position
    :param max_cov: Given max_cov from whatshap call
    :return: integer which is the maximum coverage in the whole interval, could maximal be the max_cov which was given
     as parameter in the whatshap call
    '''
    maximum_cov=max_cov
    return maximum_cov


