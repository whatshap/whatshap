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

#TODO Need to assert somewhere that if less Reads than coverage...?


def __pq_construction_out_of_given_reads( priorityqueue,readset,positions,SNP_read_map):
	'''Constructiong of the priority queue for the readset, so that each read is in the priority queue and
	sorted by their score'''

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
		for m in range(0, len(SNPset)):
			SNP_read_map[SNPset[m]].append(pq_item)

	return (priorityqueue,SNP_read_map,phasable_SNPs,vcf_indices)

def slice_read_selection(pq,coverages,already_covered_SNPS,selected_reads,phasable_SNPs,MAX_cvo,readset,Vcf_indices,SNP_read_map):
	'''Extraction of a set of read indices, where each SNP should be covered at least once, if coverage, or reads are allowing it '''
	#reads selected in this run
	#actual_selection=[]

	#TODO add an additional condition like if number covered SNPS equal phasable SNPS
	while len(pq)!= 0:
		#TODO ISEMPTY does not work

		(max_score,max_item)= pq.pop()

		if not max_item in selected_reads:
			extracted_read=readset[max_item]
			uncovered_SNP= False
			#look if positions covered by this reads are already covered or not
			for pos in extracted_read:
				#TODO could be done easier
				if pos in already_covered_SNPS:
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
						else:
							print('Current element not anymore member of the priority queue  ')

	return  (selected_reads,coverages)


def readselection(readset, positions,max_cov):
	'''The whole readselection should work in this method'''

	SNP_read_map = []
	for i in range(0, len(positions)):
		SNP_read_map.append([])
	pq=PriorityQueue()

	#Construction of SNP_read mapping and priority queue
	(pq,SNP_read_map,phasbale_SNPs,vcf_indices)=__pq_construction_out_of_given_reads( pq,readset,positions,SNP_read_map)

	#Beginning to select the reads

	#Initialize the Cov Monitor and the list of read indices which are selected..
	coverages = CovMonitor(len(positions))

	#TODO Look if we could also use there a set instead of an array
	selected_reads= []
	#Use as set
	already_covered_SNPs=set()


	(sliced_selected_reads,coverages)= slice_read_selection(pq,coverages,already_covered_SNPs,selected_reads,phasbale_SNPs,max_cov,readset,vcf_indices,SNP_read_map)
	print('selected_reads')
	print(selected_reads)
	new_set=readset
	print(readset)
	print(len(readset))
	new_read_index_set= [i for i in range(0,len(readset))]
	print('new_read_index_set')
	print(new_read_index_set)
	#Now intersection

	#(sliced_selected_reads,coverages)= slice_read_selection(pq,coverages,already_covered_SNPs,selected_reads,phasbale_SNPs,max_cov,readset,vcf_indices,SNP_read_map)









	#gives only one sclice of the reads (so where every read should be covered at least once)
	loop_integer= max_cov
	#while loop_integer!=0 :
	#	#TODO have to do :
	#	#slice the reads
	#	#do the block finding
	#	#build up a new pq with the readset without the selected reads
	#
	#
	#	(sliced_selected_reads,coverages)= slice_read_selection(pq,coverages,already_covered_SNPs,selected_reads,phasbale_SNPs,max_cov,readset,vcf_indices,SNP_read_map)
	#	loop_integer -=1


	#out= slice_read_selection(pq,coverages,already_covered_SNPs,selected_reads,phasbale_SNPs,max_cov,readset,vcf_indices,SNP_read_map)


#TODO insert the bridging .....

	#while True:
		#Need to find out if pq changes in this method or if by using this method the pq has changed, therefore we
		#need to construct the pq again (maybe without the already selected reads)

		#slice_selection(pq,coverages,already_covered_SNPs,selected_reads,phasbale_SNPs,max_cov,readset,vcf_indices,SNP_read_map)




	return selected_reads