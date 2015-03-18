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

from whatshap._core import PyRead as Read
from whatshap._core import PyReadSet as ReadSet
from whatshap._core import PyDPTable as DPTable
from whatshap._core import PyIndexSet as IndexSet
#from whatshap.scripts.whatshap import CoverageMonitor as CovMonitor
from whatshap.priorityqueue import PriorityQueue
from whatshap.coverage import CovMonitor
from .graph import ComponentFinder

#TODO Need to assert somewhere that if less Reads than coverage...?


def __pq_construction_out_of_given_reads( priorityqueue,readset,positions,SNP_read_map,bool_SNP_const):
	'''Constructiong of the priority queue for the readset, so that each read is in the priority queue and
	sorted by their score, and the given boolean says if the SNP_read_map has to be constructed.'''

	#dictionary to map the SNP positions to an index
	vcf_indices = {position: index for index, position in enumerate(positions)}


	skipped_reads = 0
	phasable_SNPs = 0
	#TODO if we want to see which SNPs are unphasable need to compute all
	#TODO  SNP positions between begin and end position.which may contribute to coverage but not to phasability

	for index, read in enumerate(readset):
		#filter out reads which cover only one variant
		if len(read) < 2:
			skipped_reads += 1
			continue
		else:
			phasable_SNPs +=1

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
		pq_item = index

		if bool_SNP_const:
			for m in range(0, len(SNPset)):
				SNP_read_map[SNPset[m]].append(pq_item)

	return (priorityqueue,SNP_read_map,phasable_SNPs,vcf_indices)

def slice_read_selection(pq,coverages,already_covered_SNPS,selected_reads,phasable_SNPs,MAX_cvo,readset,Vcf_indices,SNP_read_map):
	'''Extraction of a set of read indices, where each SNP should be covered at least once, if coverage, or reads are allowing it '''
	#Intern list for storing the actual selected reads
	intern_list=[]
	#TODO add an additional condition like if number covered SNPS equal phasable SNPS then big pq are reduced
	while len(pq)!= 0:
		#TODO ISEMPTY does not work
		(max_score,max_item)= pq.pop()
		if not max_item in selected_reads:
			extracted_read=readset[max_item]
			uncovered_SNP= False
			#look if positions covered by this reads are already covered or not
			for pos in extracted_read:
				if pos.position  in already_covered_SNPS:
					continue
				else:
					uncovered_SNP=True
					break
			#only if at least one position is not covered so the boolean is true then we could add the read if he suits into the coverage
			#Need for begin and end the vcf_index and not the correct position
			#TODO therefore maybe change the Coverage_MOnitor int the coverage.py
			begin=Vcf_indices.get(extracted_read[0].position)
			end= Vcf_indices.get(extracted_read[len(extracted_read)-1].position)
			if uncovered_SNP and coverages.max_coverage_in_range(begin,end)<MAX_cvo:
				coverages.add_read(begin,end)
				#selected_reads includes only the read indices of the selected reads  where max_item is the index....
				selected_reads.append(max_item)
				intern_list.append(max_item)

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
	return  (intern_list,selected_reads,coverages)

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


def find_new_blocks(actual_reads,all_possible_reads,component_finder,reads):
	'''Insert the new found reads into the Component finder and returns it and the corresponding components'''
	for i in actual_reads:
		read= reads[i]
		read_positions= [variant.position for variant in read if variant.position in all_possible_reads]
		for position in read_positions[1:]:
			component_finder.merge(read_positions[0], position)
	components={position : component_finder.find(position) for position in all_possible_reads}
	return (component_finder,components)


def init_blocking(vcf_indices):
	'''Initialization of the Component finder with all phasable SNPS'''
	all_possible_positions= vcf_indices.keys()
	#initializations done elsewhere because here it should be iteratively adding the new found reads into the already
	#existing component finder data structure....
	component_finder=ComponentFinder(all_possible_positions)
	return (component_finder,all_possible_positions)




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



def new_bridging (com_values,readset,selected_reads,component_finder):
	'''
	:param com_values: Components from the actual component_Finder
	:param readset: original Readset
	:param selected_reads: at the moment selected read indices
	:param component_finder: actual component finder
	:return: list of read indices which build up bridges between the given components
	'''
	outlist =[]
	#Here only looking at the positions themselves
	for index, read in enumerate(readset):
		covers_positions= 0
		#looks if read has a bridge and is not already selected
		for pos in read :
			if component_finder.find(pos.position) in com_values and index not in selected_reads:
				covers_positions+= 1

		if covers_positions>= 2:
			outlist.append(index)

	return outlist



def analyse_bridging_reads(bridging_reads, readset, selected_reads,component_finder,Cov_Monitor,vcf_indices,components,max_cov):
	''' looks at the extracted bridging reads and only select those which suit into the Coverage Monitor and only 1 read for each bridge	'''
	found_bridge= set()
	selction = []
	for index in bridging_reads:
		read= readset[index]
		for pos in read:
			read_positions=[]
			for comp in components:
				if pos.position == comp and pos.position not in found_bridge:
					found_bridge.add(pos.position)
					selction.append(index)
					begin=vcf_indices.get(read[0].position)
					end=vcf_indices.get(read[len(read)-1].position)
					if pos.position in vcf_indices.keys():
						read_positions.append(pos.position)
					if Cov_Monitor.max_coverage_in_range(begin, end ) < max_cov:
						Cov_Monitor.add_read(begin,end)
						selected_reads.append(index)
			for pos in read_positions[1:]:
				component_finder.merge(read_positions[0],pos)

	new_components={position : component_finder.find(position) for position in vcf_indices.keys()}

	return selction








def readselection(readset, positions,max_cov):
	'''The whole readselection should work in this method'''

	SNP_read_map = []
	for i in range(0, len(positions)):
		SNP_read_map.append([])



	pq=PriorityQueue()

	#Construction of SNP_read mapping and priority queue
	(pq,SNP_read_map,phasbale_SNPs,vcf_indices)=__pq_construction_out_of_given_reads( pq,readset,positions,SNP_read_map,True)

	#Initialization of Coverage Monitor
	coverages = CovMonitor(len(positions))

	#TODO Look if we could also use there a set instead of an array not yet used
	final_selected_reads=ReadSet()

	#break parameter for the loop
	#TODO maybe better parameter could be found
	loop_integer= max_cov

	#TODO change the pure index Set to a reduced read set to reduce runtime.
	#If done : 	#TODO - Has to know that the indices are shifting now, because the readset has differed....
	new_readset_subset=readset

	#TODO Look if we could also use there a set instead of an array
	selected_reads= []

	final_read_set= set()
	#Initialization of the Somponent_finder
	(component_finder,all_possible_positions)=init_blocking(vcf_indices)


	while loop_integer!=0 :
	# TODO maybe not needed to minimize the qp?
		#set includes the covered positions
		already_covered_SNPs=set()
		#TODO have to differ between REads added only in this round and all reads till there but has to be done in slice_read_selection
		(actual_reads,selected_reads,coverages)= slice_read_selection(pq,coverages,already_covered_SNPs,selected_reads,phasbale_SNPs,max_cov,new_readset_subset,vcf_indices,SNP_read_map)
		#TODO maybe NOT NEEDED because of checkout what was already selected
		#final_selected_reads=building_up_final_readset(new_readset_subset,sliced_selected_reads,final_selected_reads,original_readset)

		for slice in selected_reads:
			final_read_set.add(slice)

		#find components
		(component_finder,components)=find_new_blocks(actual_reads,all_possible_positions,component_finder,readset)
		com_values= set(components.values())

		#TODO does not work
		#New_reads=find_bridging_read(com_values,SNP_read_map,readset,actual_reads,component_finder,vcf_indices)
		bridges_reads=new_bridging(com_values,readset,selected_reads,component_finder)

		analyse_bridging_reads(bridges_reads, readset, selected_reads,component_finder,coverages,vcf_indices, SNP_read_map,actual_reads)


		#TODO SAME AS ABOVE NOT NEEDED
		#new_readset_subset = readset.subset(create_new_readset(readset,sliced_selected_reads))

		#Construction of SNP_read mapping and priority queue without SNP_map
		(pq,SNP_read_map,phasbale_SNPs,vcf_indices)=__pq_construction_out_of_given_reads( pq,new_readset_subset,positions,SNP_read_map,False)

		loop_integer -=1
		print('Components at the end')
		print(components)
		print(set(components.keys()))
		print(len(set(components.keys())))
		print(set(components.values()))
		print(len(set(components.values())))

	return selected_reads