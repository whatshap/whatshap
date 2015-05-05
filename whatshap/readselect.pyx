# THIS is totaly not used

import math
import logging
#TODO at the moment no possibility to have a cplusplus type of defaultdict
#from collections import defaultdict


from whatshap.priorityqueue import PriorityQueue
from whatshap.coverage import CovMonitor
from .graph import ComponentFinder
from whatshap.priorityqueue cimport PriorityQueue
#pyximport priorityqueue.pyx

from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.pair cimport pair
from libcpp cimport bool
#from libcpp.ordered_map cimport ordered_map
from libcpp.string cimport string
from libcpp.set cimport set


ctypedef vector[int] priority_type
ctypedef priority_type* priority_type_ptr
ctypedef int index_type
#For a read we have position,allele and quality score
#For this it does not work
#ctypedef
#ctypedef vector[int, string, int] item
#ctypedef vector[item] read_type
ctypedef vector[int] SNP_positions_type
#ctypedef vector[read_type] read_set_type
#ctypedef pair[priority_type_ptr,item_type] queue_entry_type
#For SNP_ read _MAP need ordere map
#ctypedef pair[int, priority_type ]  variant_map_type
#TODO at the moment specially for the tuple score
ctypedef vector[int] score_type
#former tuple <int, int,int> did not work

#TODO produces endless loop

logger = logging.getLogger(__name__)

#TODO Need to assert somewhere that if less Reads than coverage...?

# ====== ReadSet ======
#cdef extern from "/whatshap.priorityqueue.cpp":
#	cdef cppclass PriorityQueue:
#		bool _vector_score_lower(priority_type_ptr, priority_type_ptr) except +
#		priority_type_ptr _pyscore_to_vector(int) except +
#		int _parent(int) except +
#		int _left_child(int) except +
#		int _right_child(int) except +
#		void __dealloc__() except +
#		#ReadSet* subset(IndexSet*) except +
#		# TODO: Check why adding "except +" here doesn't compile
#		#vector[unsigned int]* get_positions()






def _construct_indexes(readset):
	''' The parameter readset: is the given ReadSet and returns, all possible variant positions
	, the vcf_index_ mapping and the SNP_read_map'''
	cdef index_type index
	cdef index_type snp_index
	cdef unordered_map[int,vector[int]] SNP_read_map
	cdef vector[int] help

	#SNP_positions_type
	positions = readset.get_positions()
	#unordered_map[int,int]
	vcf_indices = {position: index for index, position in enumerate(positions)}
	for index, read in enumerate(readset):
		for variant in read:
			snp_index = vcf_indices[variant.position]
			help= SNP_read_map[snp_index]
			help.append(index)
			SNP_read_map[snp_index] = help
			#SNP_read_map[snp_index].addindex)

	return positions, vcf_indices, SNP_read_map


def _update_score_for_reads(score_type former_score, readset, int index, already_covered_SNPs):
	'''updatest the score of the read, depending on how many reads are already covered'''
	cdef int pos
	cdef int first_score
	cdef int second_score
	cdef int quality
	#cdef read_type read
	(first_score, second_score, quality) = former_score
	read = readset[index]
	for pos in read:
		if pos.position not in already_covered_SNPs:
			first_score -= 1
	return (first_score, second_score, quality)


cdef score_type _compute_score_for_read(readset,index_type  index, vcf_indices):
	'''Method for computing the score for a read independently'''
#	#TODO At the moment tuple (new- bad ,good -bad, min(qualities))
#	cdef read_type read
	cdef int min_quality
	cdef int good_score
	cdef int bad_score
#	cdef SNP_positions_type covered_SNPS
	read = readset[index]
	#TODO look if the initial values is reasonable
	#initialize minimal quality by high value, good_score by 0  and bad_score by 0 .
	min_quality = 1000
	good_score = 0
	bad_score = 0
	covered_SNPS = []
	for pos in read:
		quality=pos.quality
		min_quality = min(min_quality, quality)
		SNP_covered = vcf_indices.get(pos.position)
		if SNP_covered != None :
			covered_SNPS.append(SNP_covered)
			#good_score represents how many SNPs are really covered
			good_score += 1
	#bad_score here means number of SNPs which are not covered...
	if len(covered_SNPS) != (covered_SNPS[len(covered_SNPS) - 1] - covered_SNPS[0] + 1):
		bad_score = (covered_SNPS[len(covered_SNPS) - 1] - covered_SNPS[0] + 1) - len(covered_SNPS)

	#initially new_score is the same as good_score, but new_score is updated later
	return (good_score - bad_score, good_score - bad_score, min_quality)


cdef priority_type __pq_construction_out_of_given_reads(readset,read_indices, vcf_indices):
	'''Constructiong of the priority queue out of the readset and the undicided read_indices,
	 so that each read is in the priority queue and sorted by their score.'''

	#Maybe combine the method with the compute score method, but for the moment the score computation should be independent
	cdef priority_type priorityqueue
	cdef int index
#	cdef read_type read
	cdef score_type computed_score
	#priorityqueue = PriorityQueue()

	for index in read_indices:
		read = readset[index]
		computed_score = _compute_score_for_read(readset, index, vcf_indices)
		#print('Computed score')

		priorityqueue.push(computed_score, index)
		#PriorityQueue.push(priorityqueue)
		#help2= (computed_score,index)
		#priorityqueue.push_back(help2)
	print('length of priorityqueue')
	print(len(priorityqueue))

	return priorityqueue


#TODO Here the pq already includes only the undecided reads....

def slice_read_selection( pq, coverages, int max_cov, readset, vcf_indices, SNP_read_map):
	'''Extraction of a set of read indices, where each SNP should be covered at least once, if coverage, or reads are allowing it '''
	#Intern list for storing the actual selected reads
#	cdef SNP_positions_type already_covered_SNPs
#	cdef set<int> reads_in_slice , reads_violating_coverage, SNPS_Covered_for_this_read, reads_whose_score_has_to_be_updated, decrease_set,d_set
	cdef int max_score
	cdef int max_item
	cdef int begin
	cdef int end
	cdef int element,
#	cdef read_type extracted_read
	cdef bool covers_new_snp
	#cdef score_type oldscore
	cdef score_type newscore
	cdef set[int] already_covered_SNPs
	cdef set[int] reads_in_slice
	cdef set[int] reads_violating_coverage
	cdef set[int] SNPS_Covered_for_this_read
	cdef set[int] reads_whose_score_has_to_be_updated

	#already_covered_SNPs =set()
	#reads_in_slice = set()
	#reads_violating_coverage = set()
	#TODO is empty does not work here...
	while len(pq) != 0 :
	#while not PriorityQueue.is_empty(pq):
		#SNPS_Covered_for_this_read=set()
		(max_score, max_item) = pq.pop()
		extracted_read = readset[max_item]
		covers_new_snp = False
		#look if positions covered by this reads are already covered or not
		for pos in extracted_read:
			it=already_covered_SNPs.find(pos.position)
			if it!= already_covered_SNPs.end():
				continue
			else:
			#if pos.position in already_covered_SNPs:
			#	continue
			#else:
				covers_new_snp = True
				#stores the positions the read covers
				SNPS_Covered_for_this_read.add(pos.position)
		#only if at least one position is not covered then we could add the read if he does not break the max coverage
		begin = vcf_indices.get(extracted_read[0].position)
		end = vcf_indices.get(extracted_read[len(extracted_read) - 1].position) + 1
		if coverages.max_coverage_in_range(begin, end) >= max_cov:
			reads_violating_coverage.add(max_item)
		elif covers_new_snp:
			coverages.add_read(begin, end)
			reads_in_slice.add(max_item)

			#TODO NOT allowed here a cdef
			#TODO aber muss hier lee erstellt werden
			#cdef set[int] reads_whose_score_has_to_be_updated
			#reads_whose_score_has_to_be_updated = set()

			#again go over the positions in the read and add them to the already_covered_SNP list

			#Only the positions in the read which cover new SNPs are analysed
			for pos in SNPS_Covered_for_this_read:
				already_covered_SNPs.add(pos)
				reads_whose_score_has_to_be_updated.update(SNP_read_map[vcf_indices.get(pos)])

			#find difference between to_decrease_score and selected_reads in order to not to try to decrease score by selected reads

			#TODO Here reads_in_slice is already a set

			#selected_read_set = set(reads_in_slice)
			decrease_set = reads_whose_score_has_to_be_updated

			#d_set = decrease_set.difference(selected_read_set)
			d_set = decrease_set.difference(reads_in_slice)

			#Catch with None if element is not in the priorityqueue
			for element in d_set:
				oldscore = pq.get_score_by_item(element)
				if oldscore != None:
					newscore = _update_score_for_reads(oldscore, readset, element, SNPS_Covered_for_this_read)
					pq.change_score(element, newscore)
	return reads_in_slice, reads_violating_coverage



def format_read_source_stats( readset, indices):
	"""Creates a string giving information on the source_ids of the reads with the given indices."""
	cdef int i
	cdef unordered_map[int,int] source_id_counts
	if len(indices) == 0:
		return 'n/a'
	#source_id_counts = defaultdict(int)
	for i in indices:
		source_id_counts[readset[i].source_id] += 1
	present_ids = list(source_id_counts.keys())
	present_ids.sort()
	return ', '.join('{}:{}'.format(source_id, count) for source_id, count in source_id_counts.items())


def readselection(readset,int  max_cov,bool bridging=True):
	'''Return the selected readindices which do not violate the maximal coverage, and additionally usage of a boolean for deciding if
	the bridging is needed or not.'''
#	cdef SNP_positions_type positions
#   cdef unordered_map 	vcf_indices
#	cdef dict SNP_read_map
	cdef set[int] selected_reads
	cdef set[int] bridging_reads
	cdef int number_unuseful_reads
	cdef int loop
	cdef int read_index
	cdef int score
	cdef int begin
	cdef int end
	cdef priority_type pq
#	cdef SNP_positions_type undecided_reads
#	cdef read_type read
#
#
#
#
	positions, vcf_indices, SNP_read_map = _construct_indexes(readset)

	logger.info('Running read selection for %d reads covering %d variants (bridging %s)', len(readset), len(positions),
				'ON' if bridging else 'OFF')

	#initialization of Coverage Monitor
	coverages = CovMonitor(len(positions))

	# indices of reads that have been selected
	#selected_reads = set()

	# indices of reads that could (potentially) still be selected ,do not consider thes read which only cover one SNP
	undecided_reads = {i for i, read in enumerate(readset) if len(read) >= 2}

	number_unuseful_reads = len(readset) - len(undecided_reads)

	loop = 0

	while loop < 15:
		print('in loop')
#	while len(undecided_reads) > 0:
		pq = __pq_construction_out_of_given_reads(readset, undecided_reads, vcf_indices)
		reads_in_slice, reads_violating_coverage = slice_read_selection(pq, coverages, max_cov, readset, vcf_indices,
																		SNP_read_map)
#		selected_reads.update(reads_in_slice)
#		undecided_reads -= reads_in_slice
#		undecided_reads -= reads_violating_coverage

		# Create new component finder from reads just selected
#		component_finder = ComponentFinder(positions)
#		for read_index in reads_in_slice:
#			read = readset[read_index]
#			read_positions = [variant.position for variant in read]
#			for position in read_positions[1:]:
#				component_finder.merge(read_positions[0], position)

		#bridging_reads = set()
#		if bridging:
#			pq = __pq_construction_out_of_given_reads(readset, undecided_reads, vcf_indices)
#			while len(pq) != 0 :
#			#TODO is_empty does not work
#			#while not pq.is_empty():
#				score, read_index = pq.pop()
#				read = readset[read_index]
#				covered_blocks = {component_finder.find(pos.position) for pos in read}

				# TODO: check coverage (potentially remove from undecided_reads)

				#Coverage Monitor
#				begin = vcf_indices.get(read[0].position)
#				end = vcf_indices.get(read[len(read) - 1].position) + 1
#				if coverages.max_coverage_in_range(begin, end) >= max_cov:
#					undecided_reads.remove(read_index)
#					continue

#				# skip read if it only covers one block
#				if len(covered_blocks) < 2:
#					continue
#				bridging_reads.add(read_index)
#				selected_reads.add(read_index)

#				coverages.add_read(begin, end)

#				undecided_reads.remove(read_index)
#				#Update Component_finder
#				read_pos = [variant.position for variant in read]
#				for pos_in_read in read_pos[1:]:
#					component_finder.merge(read_pos[0], pos_in_read)
		loop += 1


#		logger.info(
#			'... iteration %d: selected %d reads (source: %s) to cover positions and %d reads (source: %s) for bridging; %d reads left undecided',
#			loop, len(reads_in_slice), format_read_source_stats(readset, reads_in_slice), len(bridging_reads),
#			format_read_source_stats(readset, bridging_reads), len(undecided_reads)
#		)

	#for Debugging
#	new_components = {position: component_finder.find(position) for position in vcf_indices.keys()}
#	stats = (number_unuseful_reads)

	return selected_reads
#	return selected_reads, new_components, stats
