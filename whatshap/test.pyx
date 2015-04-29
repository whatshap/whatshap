from priorityqueue import PriorityQueue
from whatshap.coverage import CovMonitor

from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.pair cimport pair
from libcpp cimport bool
from libcpp cimport set

ctypedef vector[int] priority_type
ctypedef priority_type* priority_type_ptr
ctypedef int item_type
ctypedef pair[priority_type_ptr,item_type] queue_entry_type






#for using method in whatshap using def and not cdef
def addition (int x,int y):
    cdef int z
    z= x+y
    pq=PriorityQueue()
    pq.push(5,z)
    #initialization of Coverage Monitor
    coverages = CovMonitor(10)
    return z


def _construct_indexes(readset):
    ''' The parameter readset: is the given ReadSet and returns, all possible variant positions, the vcf_index_ mapping
    and the SNP_read_map'''
    cdef unordered_map[int,vector[int]] SNP_read_map
    positions = readset.get_positions()
    vcf_indices = {position: index for index, position in enumerate(positions)}

    #SNP read map with defaultdict does not work
    #SNP_read_map = defaultdict(list)
    for index, read in enumerate(readset):
        for variant in read:
            snp_index = vcf_indices[variant.position]
            #Try to avoid the append because this does not work
            diff_vector= SNP_read_map[snp_index]
            diff_vector.push_back(index)
            SNP_read_map[snp_index]= diff_vector
            #SNP_read_map[snp_index].append(index)
    return positions, vcf_indices, SNP_read_map




def readselection_2(readset, int max_cov, bridging=True):
    '''Return the selected readindices which do not violate the maximal coverage, and additionally usage of a boolean for deciding if
     the bridging is needed or not.'''
    cdef int number_unuseful_reads, lala
    #That does not work
    #cdef set[int] selected_reads

    positions, vcf_indices, SNP_read_map = _construct_indexes(readset)

	#logger.info('Running read selection for %d reads covering %d variants (bridging %s)', len(readset), len(positions),
	#'ON' if bridging else 'OFF')

	#initialization of Coverage Monitor
	#coverages = CovMonitor(len(positions))

	# indices of reads that have been selected
    selected_reads = set()

	# indices of reads that could (potentially) still be selected ,do not consider thes read which only cover one SNP
    undecided_reads = set(i for i, read in enumerate(readset) if len(read) >= 2)

    number_unuseful_reads = len(readset) - len(undecided_reads)
    print('number unuseful reads')
    print(number_unuseful_reads)

    lala=10
    return lala



