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

		#score for the sorting of best reads
		score = 0
		#set containing for each read which variants are covered by the read
		SNPset = []

		#look in the dictionary if the found position corresponds to a variantif yes increase the score
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
	#reads selected in this run
	#actual_selection=[]
	#Intern list for storing the actual selected reads
	intern_list=[]
	#TODO add an additional condition like if number covered SNPS equal phasable SNPS
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
				#selected_reads includes only the read indices of the selected reads ....
				selected_reads.append(max_item)
				intern_list.append(max_item)

				#again go over the positions in the read and add them to the already_covered_SNP list
				for pos in extracted_read:
					already_covered_SNPS.add(pos.position)

					#for extracted read decrease score of every other read which covers one of the other SNPS.
					to_decrease_score=SNP_read_map[Vcf_indices.get(pos.position)]

					#find difference between to_decrease_score and  selected_reads
					#changing to set
					#TODO maybe also possible to do this in the first place
					selected_read_set = set(selected_reads)
					decrease_set = set(to_decrease_score)

					#NEED EXACTLY TO HAVE THE ELEMENT OF THE SNP_READ_MAP WIthoiut the already selected reads
					d_set=decrease_set.difference(selected_read_set)

					#TODO need to look if element is in the heap at all ...... Catched it with None
					for element in d_set:
						oldscore=pq.get_score_by_item(element)
						if oldscore != None:
							pq.change_score(element,oldscore-1)
	return  (intern_list,selected_reads,coverages)


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
		#same code as in find_components....
		#need here a list and not a se with {}
		read_positions= [variant.position for variant in read if variant.position in all_possible_reads]
		#print('Read positions')
		#print(read_positions)
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

	#print('component_finder')
	#print(component_finder)
	#print('VCF positions ')
	#print(vcf_indices.keys())
	#for i in actual_reads:
	#	read= reads[i]
		#same code as in find_components....
		#need here a list and not a se with {}
	#	read_positions= [variant.position for variant in read if variant.position in all_possible_positions]
	#	print('Read positions')
	#	print(read_positions)
	#	for position in read_positions[1:]:
	#		component_finder.merge(read_positions[0], position)

		#print(reads[i])
	#components={position : component_finder.find(position) for position in all_possible_positions}
	#print('Components')
	#print(components)

	#TODO it does not work here because a component _finder is always new initialized....

	return (component_finder,all_possible_positions)









#How the componentfind works....
'''
def find_components(superreads, reads):
	"""
	Return a dict that maps each position to the component it is in. A
	component is identified by the position of its leftmost variant.
	"""
	logger.debug('Finding connected components ...')

	# The variant.allele attribute can be either 0 (major allele), 1 (minor allele),
	# or 3 (equal scores). If all_heterozygous is on (default), we can get
	# the combinations 0/1, 1/0 and 3/3 (the latter means: unphased).
	# If all_heterozygous is off, we can also get all other combinations.
	# In both cases, we are interested only in 0/1 and 1/0.
	phased_positions = [ v1.position for v1, v2 in zip(*superreads)
		if (v1.allele, v2.allele) in ((0, 1), (1, 0))
	]
	assert phased_positions == sorted(phased_positions)

	# Find connected components.
	# A component is identified by the position of its leftmost variant.
	component_finder = ComponentFinder(phased_positions)
	phased_positions = set(phased_positions)
	for read in reads:
		positions = [ variant.position for variant in read if variant.position in phased_positions ]
		for position in positions[1:]:
			component_finder.merge(positions[0], position)
	components = { position : component_finder.find(position) for position in phased_positions }
	logger.info('No. of variants considered for phasing: %d', len(superreads[0]))
	logger.info('No. of variants that were phased: %d', len(phased_positions))
	return components

def best_case_blocks(reads):
	"""
	Given a list of core reads, determine the number of phased blocks that
	would result if each variant were actually phased.

	Return the number of connected components.
	"""
	positions = set()
	for read in reads:
		for variant in read:
			positions.add(variant.position)
	component_finder = ComponentFinder(positions)
	for read in reads:
		read_positions = [ variant.position for variant in read ]
		for position in read_positions[1:]:
			component_finder.merge(read_positions[0], position)
	# A dict that maps each position to the component it is in.
	components = { component_finder.find(position) for position in positions }
	return len(components)




'''



















def find_bridging_read(com_values,SNP_read_map,readset,actual_reads,component_finder,vcf_indices):
	#Testing if double occurence work
	Test_set1_= set([1,2,3,4,5])
	Test_set2_=set([6,7,8,9,4])
	Check_set= Test_set1_ & Test_set2_
	#print('Test')
	#print(Check_set)
	#print('SNP Read Map')

	#print(SNP_read_map)
	all_list= [item for sublist in SNP_read_map for item in sublist]
	print('all_list.sort()')
	print(all_list)

	newlist=[]
	double_list=[]
	for i in all_list:
		if i not in newlist:
			newlist.append(i)
		else:
			double_list.append(i)
	print('Double List')
	print(double_list)
	print(len(double_list))
	print(len(set(double_list)))
	print('Newlist')
	print(newlist)
	print(len(newlist))

	#TODO At the moment in both sets the double and the newlist the length is equal this means there is no read which covers
	#TODO 2 representatives and is not added already to the components


	new_reads= []
	#print('FIND BRIDGING')
	representatives=[vcf_indices[val] for val in com_values]
	#print('Representatives')
	#print(representatives)
	#alle möglichen reads rausfinden die die Representatnen abdecken, dabei die ohne Reads rauswerfen
	SNP_map_of_reps= [SNP_read_map[vcf] for vcf in representatives if SNP_read_map[vcf]!= []]

	#print('All lists')
	#print(all_list)

	#print('SNP_map_of_reps')
	#print(SNP_map_of_reps)

	Double_occurence_set= set()

#TODO Could not work in this because one set is empty
	for reas_list in SNP_map_of_reps:
		new_set=set(reas_list)
		Double_occurence_set= Double_occurence_set & new_set

		#print('Double occurence set')
		#print(Double_occurence_set)


		#print(val)
		#print(vcf_indices[val])
		#print(SNP_read_map[vcf_indices[val]])
		#print(SNP_read_map)
		lal=0
	return new_reads











def readselection(readset, positions,max_cov):
	'''The whole readselection should work in this method'''

	SNP_read_map = []
	for i in range(0, len(positions)):
		SNP_read_map.append([])



	pq=PriorityQueue()

	#Construction of SNP_read mapping and priority queue
	(pq,SNP_read_map,phasbale_SNPs,vcf_indices)=__pq_construction_out_of_given_reads( pq,readset,positions,SNP_read_map,True)

	#Beginning to select the reads

	#Initialize the Cov Monitor and the list of read indices which are selected..
	coverages = CovMonitor(len(positions))

	#TODO Look if we could also use there a set instead of an array

	final_selected_reads=ReadSet()

	#gives only one sclice of the reads (so where every read should be covered at least once)
	loop_integer= max_cov

	new_readset_subset=readset
	#TODO Look if we could also use there a set instead of an array
	selected_reads= []

	final_read_set= set()



	#Loop where the slicing read selection and the corresponding bridging is done


	#TODO - Has to know that the indices are shifting now, because the readset has differed....


	#TODO Do the initialization before the loop
	(component_finder,all_possible_positions)=init_blocking(vcf_indices)
	print('In while loop SNP READ MAP AFTER INITIALIZATION')
	print(SNP_read_map)

	#For testing  only one iteration
	#loop_integer=1

	while loop_integer!=0 :
		print('In while loop SNP READ MAP AFTER INITIALIZATION')
		print(SNP_read_map)



	#	#TODO have to do :
	#	#slice the reads
	#	#do the block finding
	#	#build up a new pq with the readset without the selected reads TODO maybe not needed to minimize the qp?
	#

		# #should always be an empty set, because for each slice it should start new

		already_covered_SNPs=set()
		former_selected_reads= selected_reads
		#TODO have to differ between REads added only in this round and all reads till there but has to be done in slice_read_selection
		(actual_reads,selected_reads,coverages)= slice_read_selection(pq,coverages,already_covered_SNPs,selected_reads,phasbale_SNPs,max_cov,new_readset_subset,vcf_indices,SNP_read_map)
		#TODO maybe not needed because of checkout what was already selected
		#final_selected_reads=building_up_final_readset(new_readset_subset,sliced_selected_reads,final_selected_reads,original_readset)

		#print('actual reads')
		#print(actual_reads)

		for slice in selected_reads:
			final_read_set.add(slice)
		#print('SELECTED READS')
		#print(selected_reads)
		#TODO here do the bridging
		(component_finder,components)=find_new_blocks(actual_reads,all_possible_positions,component_finder,readset)

		print('Components')
		print(components.keys())
		com_values= set(components.values())
		print(com_values)
		#for val in com_values:
		#print(val)
		#print(vcf_indices[val])
		#print(SNP_read_map[vcf_indices[val]])
			#print(SNP_read_map)
		New_reads=find_bridging_read(com_values,SNP_read_map,readset,actual_reads,component_finder,vcf_indices)
		#TODO SAME AS ABOVE NOT NEEDED
		#new_readset_subset = readset.subset(create_new_readset(readset,sliced_selected_reads))





		#final_selected_reads.add(tuple(sliced_selected_reads))
		#Construction of SNP_read mapping and priority queue
		(pq,SNP_read_map,phasbale_SNPs,vcf_indices)=__pq_construction_out_of_given_reads( pq,new_readset_subset,positions,SNP_read_map,False)

		#print('Coverage')
		#print(coverages.max_coverage_in_range(0,len(vcf_indices)))
		loop_integer -=1






	#print('complete read set')
	#print(final_read_set)
	#print('Length of the complete Read SEt')
	#print(len(selected_reads))

	#out= slice_read_selection(pq,coverages,already_covered_SNPs,selected_reads,phasbale_SNPs,max_cov,readset,vcf_indices,SNP_read_map)


	return selected_reads