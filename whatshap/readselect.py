# TODO
# Class Reads where access through SNP position -DONE in SNP MAP
# implement readscore - DONE in combined method to build up SNP MAP which include the readscore
# Heap implemented for storage of Reads - Done

#include the Coverage Monitor...

#Looking if heapq is possible to manage the heap or not especially  id the runtime changes,,
# Redefine ComponentFinder by using max value and also the rank (not sure)
#implement heuristic
#return Readset in the former representation
#Erase SNPS which are unphasable -Not Erase out of structure but counted

import math
import logging
from collections import defaultdict

from whatshap._core import PyRead as Read
from whatshap._core import PyReadSet as ReadSet
from whatshap._core import PyDPTable as DPTable
from whatshap._core import PyIndexSet as IndexSet
#from whatshap.scripts.whatshap import CoverageMonitor as CovMonitor
from whatshap.priorityqueue import PriorityQueue
from whatshap.coverage import CovMonitor
from .graph import ComponentFinder

logger = logging.getLogger(__name__)

#TODO Need to assert somewhere that if less Reads than coverage...?

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
			
		
def __pq_construction_out_of_given_reads(readset, read_indices, vcf_indices):
	'''Constructiong of the priority queue out of the readset and the undicided read_indices,
	 so that each read is in the priority queue and sorted by their score.'''

	priorityqueue = PriorityQueue()

	#TODO  SNP positions between begin and end position.which may contribute to coverage but not to phasability

	for index in read_indices:
		read = readset[index]

		#score for the sorting of best reads
		score = 0
		#set containing for each read which variants are covered by the read
		SNPset = []

		#look in the dictionary if the found position corresponds to a variant if yes increase the score
		for pos in read:
			variant_index=vcf_indices.get(pos.position)
			if variant_index!= None:
				SNPset.append(variant_index)
				score += 1

		#decrease score if SNPs covered physically, but are not sequenced (e.g. paired_end)
		if len(SNPset)!= (SNPset[len(SNPset)-1] - SNPset[0]+1):
			score = score - ((SNPset[len(SNPset)-1] - SNPset[0]+1)-len(SNPset))

		priorityqueue.push(score,index)

	return priorityqueue

#TODO Here the pq already includes only the undecided reads....

def slice_read_selection(pq,coverages,max_cov,readset,vcf_indices,SNP_read_map):
	'''Extraction of a set of read indices, where each SNP should be covered at least once, if coverage, or reads are allowing it '''
	#Intern list for storing the actual selected reads
	already_covered_SNPs=set()
	# reads that are seleted in this slice
	reads_in_slice = set()
	reads_violating_coverage = set()
	#TODO add an additional condition like if number covered SNPS equal phasable SNPS then big pq are reduced unnecessary because pq already includes only the undicided ones....
	while not pq.is_empty():
		(max_score,max_item)= pq.pop()
		extracted_read=readset[max_item]
		covers_new_snp = False
		#look if positions covered by this reads are already covered or not
		for pos in extracted_read:
			if pos.position  in already_covered_SNPs:
				continue
			else:
				covers_new_snp=True
				break
		#only if at least one position is not covered then we could add the read if he does not break the max coverage
		begin=vcf_indices.get(extracted_read[0].position)
		end= vcf_indices.get(extracted_read[len(extracted_read)-1].position) +1
		if coverages.max_coverage_in_range(begin,end) >= max_cov:
			reads_violating_coverage.add(max_item)
		elif covers_new_snp:
			coverages.add_read(begin,end)
			reads_in_slice.add(max_item)

			#again go over the positions in the read and add them to the already_covered_SNP list
			for pos in extracted_read:
				already_covered_SNPs.add(pos.position)

				#for extracted read decrease score of every other read which covers one of the other SNPS.
				to_decrease_score=SNP_read_map[vcf_indices.get(pos.position)]

				#find difference between to_decrease_score and selected_reads in order to not to try to decrease score by selected reads
				selected_read_set = set(reads_in_slice)
				decrease_set = set(to_decrease_score)
				d_set=decrease_set.difference(selected_read_set)

				#Catch with None if element is not in the priorityqueue
				for element in d_set:
					oldscore=pq.get_score_by_item(element)
					if oldscore != None:
						pq.change_score(element,oldscore-1)
	return reads_in_slice, reads_violating_coverage

#Not Needed
def new_bridging (readset,selected_reads,component_finder):
	'''
	:param readset: original Readset
	:param selected_reads: at the moment selected read indices
	:param component_finder: actual component finder
	:return: list of read indices which build up bridges between the given components
	'''
	outlist = []
	#Here only looking at the positions themselves
	for index, read in enumerate(readset):
		covered_blocks = set(component_finder.find(pos.position) for pos in read)
		
		# skip read if it only covers one block
		if len(covered_blocks) < 2:
			continue
		
		# skip read if it has already been selected
		if index in selected_reads:
			continue

		outlist.append(index)

	return outlist


#Not needed
def analyse_bridging_reads(bridging_reads, readset, selected_reads,component_finder,Cov_Monitor,vcf_indices,components,max_cov):
	'''	looks at the extracted bridging reads and only select those which suit into the Coverage Monitor and only 1 read for each bridge	'''

	found_bridge= set()
	selction = []
	for index in bridging_reads:
		read= readset[index]
		for pos in read:
			read_positions=[]
			for comp in components:
				#TODO Here Maybe find more than one bridge to one component .....
				if pos.position == comp and pos.position not in found_bridge:
					found_bridge.add(pos.position)
					selction.append(index)
					begin=vcf_indices.get(read[0].position)
					end=vcf_indices.get(read[len(read)-1].position)+1
					if pos.position in vcf_indices.keys():
						read_positions.append(pos.position)
					if Cov_Monitor.max_coverage_in_range(begin, end ) < max_cov:
						#print('Added read index')
						#print(index)
						Cov_Monitor.add_read(begin,end)

						selected_reads.add(index)
			for pos in read_positions[1:]:
				component_finder.merge(read_positions[0],pos)

	new_components={position : component_finder.find(position) for position in vcf_indices.keys()}

	#print('New Components after bridging')
	#print(set(new_components.values()))
	#print(len(set(new_components.values())))

	return (selction,Cov_Monitor,selected_reads,component_finder)


def format_read_source_stats(readset, indices):
	"""Creates a string giving information on the source_ids of the reads with the given indices."""
	if len(indices) == 0:
		return 'n/a'
	source_id_counts = defaultdict(int)
	for i in indices:
		source_id_counts[readset[i].source_id] += 1
	present_ids = list(source_id_counts.keys())
	present_ids.sort()
	return ', '.join('{}:{}'.format(source_id,count) for source_id,count in source_id_counts.items())


def readselection(readset, max_cov, bridging = True):
	'''Return the selected readindices which do not violate the maximal coverage, and additionally usage of a boolean for deciding if
	 the bridging is needed or not.'''

	positions, vcf_indices, SNP_read_map = _construct_indexes(readset)

	logger.info('Running read selection for %d reads covering %d variants (bridging %s)', len(readset), len(positions), 'ON' if bridging else 'OFF')

	#initialization of Coverage Monitor
	coverages = CovMonitor(len(positions))

	# indices of reads that have been selected
	selected_reads = set()

	# indices of reads that could (potentially) still be selected ,do not consider thes read which only cover one SNP
	undecided_reads = set( i for i, read in enumerate(readset) if len(read) >= 2 )
	number_unuseful_reads = len(readset)- len(undecided_reads)


	loop = 0
	while len(undecided_reads) > 0:
		pq = __pq_construction_out_of_given_reads(readset, undecided_reads, vcf_indices)
		reads_in_slice, reads_violating_coverage = slice_read_selection(pq,coverages,max_cov,readset,vcf_indices,SNP_read_map)
		selected_reads.update(reads_in_slice)
		undecided_reads -= reads_in_slice
		undecided_reads -= reads_violating_coverage

		# Create new component finder from reads just selected
		component_finder = ComponentFinder(positions)
		for read_index in reads_in_slice:
			read = readset[read_index]
			read_positions = [variant.position for variant in read]
			for position in read_positions[1:]:
				component_finder.merge(read_positions[0], position)

		bridging_reads = set()
		if bridging:
			pq = __pq_construction_out_of_given_reads(readset, undecided_reads, vcf_indices)
			while not pq.is_empty():
				score, read_index = pq.pop()
				read = readset[read_index]
				covered_blocks = set(component_finder.find(pos.position) for pos in read)
				
				# TODO: check coverage (potentially remove from undecided_reads)

				#Coverage Monitor
				begin=vcf_indices.get(read[0].position)
				end=vcf_indices.get(read[len(read)-1].position)+1
				if coverages.max_coverage_in_range(begin, end ) >= max_cov:
					undecided_reads.remove(read_index)
					continue

				# skip read if it only covers one block
				if len(covered_blocks) < 2:
					continue
				bridging_reads.add(read_index)
				selected_reads.add(read_index)

				coverages.add_read(begin,end)

				undecided_reads.remove(read_index)
				#Update Component_finder
				read_pos= [variant.position for variant in read]
				for pos_in_read in read_pos[1:]:
					component_finder.merge(read_pos[0],pos_in_read)
		loop +=1
		logger.info('... iteration %d: selected %d reads (source: %s) to cover positions and %d reads (source: %s) for bridging; %d reads left undecided', 
			loop, len(reads_in_slice), format_read_source_stats(readset, reads_in_slice), len(bridging_reads), format_read_source_stats(readset, bridging_reads), len(undecided_reads)
		)

	#for Debugging
	new_components={position : component_finder.find(position) for position in vcf_indices.keys()}
	#print('New componentes')
	#print(new_components)
	#print('Length of components')
	#print(len(new_components))

	#statistic for logger in whatshap:
    #1. number_unuseful_reads do skipped reads
	stats = (number_unuseful_reads)

	return selected_reads, new_components, stats
