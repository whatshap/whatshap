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
from collections import defaultdict

from whatshap._core import PyRead as Read
from whatshap._core import PyReadSet as ReadSet
from whatshap._core import PyDPTable as DPTable
from whatshap._core import PyIndexSet as IndexSet
#from whatshap.scripts.whatshap import CoverageMonitor as CovMonitor
from whatshap.priorityqueue import PriorityQueue
from whatshap.coverage import CovMonitor
from .graph import ComponentFinder

#TODO Need to assert somewhere that if less Reads than coverage...?

def _construct_indexes(readset):
	positions = readset.get_positions()
	vcf_indices = {position: index for index, position in enumerate(positions)}
	SNP_read_map = defaultdict(list)
	for index, read in enumerate(readset):
		for variant in read:
			snp_index = vcf_indices[variant.position]
			SNP_read_map[snp_index].append(index)
	return positions, vcf_indices, SNP_read_map
			
		

def __pq_construction_out_of_given_reads(readset, read_indices, vcf_indices):
	'''Constructiong of the priority queue for the readset, so that each read is in the priority queue and
	sorted by their score, and the given boolean says if the SNP_read_map has to be constructed. Will only
	use reads whose indices are given in read_indices (and will ignores the rest).'''

	priorityqueue = PriorityQueue()

	#TODO if we want to see which SNPs are unphasable need to compute all
	#TODO  SNP positions between begin and end position.which may contribute to coverage but not to phasability

	for index in read_indices:
		read = readset[index]

	#TODO decrease the readset by the not necessary reads which cover only 1 or which are already selected

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

def slice_read_selection(pq,coverages,MAX_cvo,readset,Vcf_indices,SNP_read_map):
	'''Extraction of a set of read indices, where each SNP should be covered at least once, if coverage, or reads are allowing it '''
	#Intern list for storing the actual selected reads
	already_covered_SNPs=set()
	# reads that are seleted in this slice
	reads_in_slice = set()
	reads_violating_coverage = set()
	#TODO add an additional condition like if number covered SNPS equal phasable SNPS then big pq are reduced
	while len(pq)!= 0:
		#TODO ISEMPTY does not work
		(max_score,max_item)= pq.pop()
		extracted_read=readset[max_item]
		covers_new_snp = False
		#look if positions covered by this reads are already covered or not
		for pos in extracted_read:
			if pos.position  in already_covered_SNPS:
				continue
			else:
				covers_new_snp=True
				break
		#only if at least one position is not covered so the boolean is true then we could add the read if he suits into the coverage
		#Need for begin and end the vcf_index and not the correct position
		#TODO therefore maybe change the Coverage_MOnitor int the coverage.py
		begin=Vcf_indices.get(extracted_read[0].position)
		end= Vcf_indices.get(extracted_read[len(extracted_read)-1].position)
		if coverages.max_coverage_in_range(begin,end) >= MAX_cvo:
			reads_violating_coverage.add(max_item)
		elif covers_new_snp:
			coverages.add_read(begin,end)
			#selected_reads includes only the read indices of the selected reads  where max_item is the index....
			selected_reads.add(max_item)
			reads_in_slice.append(max_item)

			#again go over the positions in the read and add them to the already_covered_SNP list
			for pos in extracted_read:
				already_covered_SNPS.add(pos.position)

				#for extracted read decrease score of every other read which covers one of the other SNPS.
				to_decrease_score=SNP_read_map[Vcf_indices.get(pos.position)]

				#find difference between to_decrease_score and selected_reads in order to not to try to decrease score by selected reads
				#TODO maybe also possible to do this in the first place by getting another readset than the origin
				selected_read_set = set(selected_reads)
				decrease_set = set(to_decrease_score)
				d_set=decrease_set.difference(selected_read_set)

				#TODO need to look if element is in the heap at all ...... Catched it with None
				for element in d_set:
					oldscore=pq.get_score_by_item(element)
					if oldscore != None:
						pq.change_score(element,oldscore-1)
	return reads_in_slice, reads_violating_coverage

#TODO at the moment not used
def create_new_readset(readset,erasing_indices):
	'''returns given a readset and a list of indices which are to erased out od the readset and returns the corresponding IndexSet'''
	new_set= IndexSet()
	readset_indices= set (i for i in range(0,len(readset)))
	for j in readset_indices:
		if j not in erasing_indices:
			new_set.add(j)
	return new_set


#TODO Setting up another readset does not work, maybe not needed
def building_up_final_readset(new_readset_subset,sliced_selected_reads,final_selected_reads,original_readset):
	'''out of a given Readset (new_readset_subset) and selected indices, build up a readset in final_selected_reads'''
	#for i in sliced_selected_reads:
	#	read= new_readset_subset[i]
	#
	#	print('Found read')
	#	print(read)
	#	final_selected_reads.push(read)
	#print('final_selected_reads')
	#print(final_selected_reads)
	return final_selected_reads




#TODO At the moment not needed
def find_bridging_read(com_values,SNP_read_map,vcf_indices):
	'''for the actual components, the readset and the componentfinder , find reads which build bridges between the components'''
	all_list= [item for sublist in SNP_read_map for item in sublist]

	newlist=[]
	double_list=[]
	for i in all_list:
		if i not in newlist:
			newlist.append(i)
		else:
			double_list.append(i)

	#TODO At the moment in both sets the double and the newlist the length is equal this means there is no read which covers
	#TODO 2 representatives and is not added already to the components


	new_reads= []
	representatives=[vcf_indices[val] for val in com_values]

	SNP_map_of_reps= [SNP_read_map[vcf] for vcf in representatives if SNP_read_map[vcf]!= []]

	Double_occurence_set= set()

	#TODO Could not work in this because one set is empty
	for reas_list in SNP_map_of_reps:
		new_set=set(reas_list)
		Double_occurence_set= Double_occurence_set & new_set
	return new_reads



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
					end=vcf_indices.get(read[len(read)-1].position)
					if pos.position in vcf_indices.keys():
						read_positions.append(pos.position)
					if Cov_Monitor.max_coverage_in_range(begin, end ) < max_cov:
						Cov_Monitor.add_read(begin,end)
						selected_reads.add(index)
			for pos in read_positions[1:]:
				component_finder.merge(read_positions[0],pos)

	new_components={position : component_finder.find(position) for position in vcf_indices.keys()}

	print('New Components after bridging')
	print(set(new_components.values()))
	print(len(set(new_components.values())))

	return (selction,Cov_Monitor,selected_reads,component_finder)








def readselection(readset, max_cov, bridging = True):
	'''The whole readselection should work in this method'''
	positions, vcf_indices, SNP_read_map = _construct_indexes(readset)

	#Initialization of Coverage Monitor
	coverages = CovMonitor(len(positions))

	#TODO Look if we could also use there a set instead of an array not yet used
	final_selected_reads=ReadSet()

	# indices of reads that have been selected
	selected_reads = set()

	# indices of reads that could (potentially) still be selected
	undecided_reads = set( i for i, read in enumerate(readset) if len(read) >= 2 )

	component_finder = ComponentFinder(positions)

	while len(undecided_reads) > 0:
		pq = __pq_construction_out_of_given_reads(readset, undecided_reads, vcf_indices)

		reads_in_slice, reads_violating_coverage = slice_read_selection(pq,coverages,max_cov,new_readset_subset,vcf_indices,SNP_read_map)
		selected_reads.update(reads_in_slice)
		undecided_reads -= reads_in_slice
		undecided_reads -= reads_violating_coverage

		# Update component finder with newly selected reads
		for read_index in reads_in_slice:
			read = readset[read_index]
			read_positions = [variant.position for variant in read]
			for position in read_positions[1:]:
				component_finder.merge(read_positions[0], position)
		#components={position : component_finder.find(position) for position in all_possible_reads}
	
		#TODO maybe NOT NEEDED because of checkout what was already selected
		#final_selected_reads=building_up_final_readset(new_readset_subset,sliced_selected_reads,final_selected_reads,original_readset)

		if bridging:
			pq = __pq_construction_out_of_given_reads(readset, undecided_reads, vcf_indices)
			while not pq.is_empty():
				score, read_index = pq.pop()
				read = readset[read_index]
				covered_blocks = set(component_finder.find(pos.position) for pos in read)
				
				# TODO: check coverage (potentially remove from undecided_reads)
				
				# skip read if it only covers one block
				if len(covered_blocks) < 2:
					continue
				
				# TODO: add to selected, remove from undecided, add to coverage monitor, update component_finder

			#print('Selection of bridging')
			#print(selction)

		#TODO SAME AS ABOVE NOT NEEDED
		#new_readset_subset = readset.subset(create_new_readset(readset,sliced_selected_reads))


		#print('Components at the end')
		#print(components)
		#print(set(components.keys()))
		#print(len(set(components.keys())))
		#print(set(components.values()))
		#print(len(set(components.values())))

	return (selected_reads,components.keys(),components.values())
