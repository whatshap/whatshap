# TODO
# Class Reads where access through SNP position -DONE in SNP MAP
# implement readscore - DONE in combined method to build up SNP MAP which include the readscore
# Heap implemented for storage of Reads - Done

# include the Coverage Monitor...
#Looking if heapq is possible to manage the heap or not especially  id the runtime changes,,
# Redefine ComponentFinder by using max value and also the rank (not sure)
#implement heuristic
#return Readset in the former representation
#Erase SNPS which are unphasable -Not Erase out of structure but counted

import math
import logging
from collections import defaultdict

from whatshap._core cimport PyRead
from whatshap._core cimport PyReadSet
#from whatshap._core import PyDPTable as DPTable
#from whatshap._core import PyIndexSet as IndexSet
#from whatshap.scripts.whatshap import CoverageMonitor as CovMonitor
from whatshap.priorityqueue cimport PriorityQueue, priority_type, priority_type_ptr, queue_entry_type
from whatshap.coverage import CovMonitor
from .graph import ComponentFinder
from libcpp.vector cimport vector
from libcpp.unordered_set cimport unordered_set

ctypedef vector[unsigned int] positions_type
ctypedef positions_type* positions_type_ptr

logger = logging.getLogger(__name__)

def _construct_indexes(readset):
	''' The parameter readset: is the given ReadSet and returns, all possible variant positions, the vcf_index_ mapping
    and the SNP_read_map'''
	positions = readset.get_positions()
	vcf_indices = {position: index for index, position in enumerate(positions)}
	SNP_read_map = defaultdict(list)
	for index, read in enumerate(readset):
		for variant in read:
			snp_index = vcf_indices[variant.position]
			SNP_read_map[snp_index].append(index)
	return positions, vcf_indices, SNP_read_map

cdef priority_type_ptr _update_score_for_reads(priority_type_ptr former_score, PyReadSet readset, index, unordered_set[int]& already_covered_SNPs):
	'''updatest the score of the read, depending on how many reads are already covered'''
	cdef int first_score = former_score.at(0)
	cdef int second_score = former_score.at(1)
	cdef int quality = former_score.at(2)
	cdef PyRead read = readset.get(index)
	cdef unordered_set[int].iterator it 
	for i in range(read.getVariantCount()):
		it = already_covered_SNPs.find(read.getPosition(i))
		if it == already_covered_SNPs.end():
			first_score -= 1
	cdef priority_type_ptr result = new priority_type()
	result.push_back(first_score)
	result.push_back(second_score)
	result.push_back(quality)
	return result


def _compute_score_depedning_on_quality_only (readset,index,vcf_indices):
	'''Method which computes anothe readscore depending only in the quality,mor precisely the
	score of the read is the average quality in the read over all SNP positions covered by the read'''
	read= readset[index]
	#min_quality= 1000
	actual_quality= 0
	for pos in enumerate:
		quality= pos.quality
		actual_quality += quality
	return (actual_quality/ len(read))



cdef priority_type_ptr _compute_score_for_read(PyReadSet readset, int index, vcf_indices):
	'''Method for computing the score for a read independently'''
	#TODO At the moment tuple (new- bad ,good -bad, min(qualities))
	cdef PyRead read = readset.get(index)
	#TODO look if the initial values is reasonable
	#initialize minimal quality by high value, good_score by 0  and bad_score by 0 .
	cdef int min_quality = 1000
	cdef int good_score = 0
	cdef int bad_score = 0
	cdef int quality = -1
	cdef int pos = -1
	covered_SNPS = []
	for i in range(read.getVariantCount()):
		quality = read.getVariantQuality(i)
		pos = read.getPosition(i)
		min_quality = min(min_quality, quality)
		SNP_covered = vcf_indices.get(pos)
		if SNP_covered != None:
			covered_SNPS.append(SNP_covered)
			#good_score represents how many SNPs are really covered
			good_score += 1
	#bad_score here means number of SNPs which are not covered...
	if len(covered_SNPS) != (covered_SNPS[len(covered_SNPS) - 1] - covered_SNPS[0] + 1):
		bad_score = (covered_SNPS[len(covered_SNPS) - 1] - covered_SNPS[0] + 1) - len(covered_SNPS)
	cdef priority_type_ptr result = new priority_type()
	#initially new_score is the same as good_score, but new_score is updated later
	result.push_back(good_score - bad_score)
	result.push_back(good_score - bad_score)
	result.push_back(min_quality)
	return result

cdef PriorityQueue __pq_construction_out_of_given_reads(PyReadSet readset, read_indices, vcf_indices):
	'''Constructiong of the priority queue out of the readset and the undicided read_indices,
     so that each read is in the priority queue and sorted by their score.'''

	#Maybe combine the method with the compute score method, but for the moment the score computation should be independent

	cdef PriorityQueue priorityqueue = PriorityQueue()

	for index in read_indices:
		read = readset.get(index)
		computed_score = _compute_score_for_read(readset, index, vcf_indices)
		priorityqueue.c_push(computed_score, index)

	return priorityqueue


#TODO Here the pq already includes only the undecided reads....

cdef slice_read_selection(PriorityQueue pq, coverages, max_cov, PyReadSet readset, vcf_indices, SNP_read_map):
	'''Extraction of a set of read indices, where each SNP should be covered at least once, if coverage, or reads are allowing it '''
	#Intern list for storing the actual selected reads
	already_covered_SNPs = set()
	reads_in_slice = set()
	reads_violating_coverage = set()
	cdef PyRead extracted_read
	cdef queue_entry_type entry
	cdef priority_type_ptr max_score
	cdef int max_item
	cdef int pos
	cdef unordered_set[int] SNPS_Covered_for_this_read 
	cdef priority_type_ptr oldscore
	cdef priority_type_ptr newscore
	while not pq.c_is_empty():
		SNPS_Covered_for_this_read.clear()
		entry = pq.c_pop()
		max_score = entry.first
		max_item = entry.second
		extracted_read = readset.get(max_item)
		covers_new_snp = False
		#look if positions covered by this reads are already covered or not
		for i in range(extracted_read.getVariantCount()):
			pos = extracted_read.getPosition(i)
			if pos in already_covered_SNPs:
				continue
			else:
				covers_new_snp = True
				#stores the positions the read covers
				SNPS_Covered_for_this_read.add(pos.position)
		#only if at least one position is not covered then we could add the read if he does not break the max coverage
		begin = vcf_indices.get(extracted_read.getPosition(0))
		end = vcf_indices.get(extracted_read.getPosition(extracted_read.getVariantCount() - 1)) + 1
		if coverages.max_coverage_in_range(begin, end) >= max_cov:
			reads_violating_coverage.add(max_item)
		elif covers_new_snp:
			coverages.add_read(begin, end)
			reads_in_slice.add(max_item)
			reads_whose_score_has_to_be_updated = set()

			#again go over the positions in the read and add them to the already_covered_SNP list

			#Only the positions in the read which cover new SNPs are analysed
			for pos in SNPS_Covered_for_this_read:
				already_covered_SNPs.add(pos)
				reads_whose_score_has_to_be_updated.update(SNP_read_map[vcf_indices.get(pos)])

			#find difference between to_decrease_score and selected_reads in order to not to try to decrease score by selected reads
			selected_read_set = set(reads_in_slice)
			decrease_set = reads_whose_score_has_to_be_updated
			d_set = decrease_set.difference(selected_read_set)

			#Catch with None if element is not in the priorityqueue
			for element in d_set:
				oldscore = pq.c_get_score_by_item(element)
				if oldscore != NULL:
					newscore = _update_score_for_reads(oldscore, readset, element, SNPS_Covered_for_this_read)
					pq.c_change_score(element, newscore)
	return reads_in_slice, reads_violating_coverage



def format_read_source_stats(readset, indices):
	"""Creates a string giving information on the source_ids of the reads with the given indices."""
	if len(indices) == 0:
		return 'n/a'
	source_id_counts = defaultdict(int)
	for i in indices:
		source_id_counts[readset[i].source_id] += 1
	present_ids = list(source_id_counts.keys())
	present_ids.sort()
	return ', '.join('{}:{}'.format(source_id, count) for source_id, count in source_id_counts.items())


def readselection(PyReadSet readset, max_cov, bridging=True):
	'''Return the selected readindices which do not violate the maximal coverage, and additionally usage of a boolean for deciding if
     the bridging is needed or not.'''
     
	positions, vcf_indices, SNP_read_map = _construct_indexes(readset)
	
	logger.info('Running read selection for %d reads covering %d variants (bridging %s)', len(readset), len(positions),'ON' if bridging else 'OFF')

	#initialization of Coverage Monitor
	coverages = CovMonitor(len(positions))

	# indices of reads that have been selected
	selected_reads = set()

	# indices of reads that could (potentially) still be selected ,do not consider thes read which only cover one SNP
	undecided_reads = set(i for i, read in enumerate(readset) if len(read) >= 2)

	number_unuseful_reads = len(readset) - len(undecided_reads)

	cdef PriorityQueue pq 
	cdef PyRead read
	cdef unordered_set[int] covered_blocks
	loop = 0
	while len(undecided_reads) > 0:
		pq = __pq_construction_out_of_given_reads(readset, undecided_reads, vcf_indices)
		reads_in_slice, reads_violating_coverage = slice_read_selection(pq, coverages, max_cov, readset, vcf_indices, SNP_read_map)
		selected_reads.update(reads_in_slice)
		undecided_reads -= reads_in_slice
		undecided_reads -= reads_violating_coverage

		# Create new component finder from reads just selected
		component_finder = ComponentFinder(positions)
		for read_index in reads_in_slice:
			read = readset.get(read_index)
			for i in range(1, read.getVariantCount()):
				component_finder.merge(read.getPosition(0), read.getPosition(i))

		bridging_reads = set()
		if bridging:
			pq = __pq_construction_out_of_given_reads(readset, undecided_reads, vcf_indices)
			while not pq.is_empty():
				score, read_index = pq.pop()
				read = readset.get(read_index)
				covered_blocks.clear()
				for i in range(read.getVariantCount()):
					covered_blocks.add(component_finder.find(read.getPosition(i)))

				# TODO: check coverage (potentially remove from undecided_reads)

				#Coverage Monitor
				begin = vcf_indices.get(read.getPosition(0))
				end = vcf_indices.get(read.getPosition(read.getVariantCount() - 1)) + 1
				if coverages.max_coverage_in_range(begin, end) >= max_cov:
					undecided_reads.remove(read_index)
					continue

				# skip read if it only covers one block
				if covered_blocks.size() < 2:
					continue
				bridging_reads.add(read_index)
				selected_reads.add(read_index)

				coverages.add_read(begin, end)

				undecided_reads.remove(read_index)
				#Update Component_finder
				for i in range(1, read.getVariantCount()):
					component_finder.merge(read.getPosition(0), read.getPosition(i))
		loop += 1
		logger.info(
			'... iteration %d: selected %d reads (source: %s) to cover positions and %d reads (source: %s) for bridging; %d reads left undecided',
			loop, len(reads_in_slice), format_read_source_stats(readset, reads_in_slice), len(bridging_reads),
			format_read_source_stats(readset, bridging_reads), len(undecided_reads)
		)

	#for Debugging
	new_components = {position: component_finder.find(position) for position in vcf_indices.keys()}
	stats = (number_unuseful_reads)

	return selected_reads, new_components, stats
