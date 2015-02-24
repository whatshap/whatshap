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


class Bestreads:
	def __init__(self, readset, positions):
		'''Initialize with the readset containing all reads and the positions which are the SNP variants'''
		self.readset = readset
		self.positions = positions


	def priorityqueue_construction(self, priorityqueue):
		'''Constructiong of the priority queue for the readset, so that each read is in the priority queue and
		sorted by their score'''
		# TODO: use list instead of dict
		SNP_read_map = {}
		for i in range(0, len(self.positions)):
			SNP_read_map[i] = []

		#dictionary to map the SNP positions to an index
		vcf_indices = {position: index for index, position in enumerate(self.positions)}


		skipped_reads = 0
		#if we want to see which SNPs are unphasable need to compute all
		# SNP positions between begin and end position.which may contribute to coverage but not to phasability

		for index, read in self.readset:


			if len(read) < 2:
				skipped_reads += 1
				continue

			#score for the sorting of best reads
			score = 0
			#set containing for each read which variants are covered by the read
			SNPset = []

			#for the reads and the vcf indices
			for j in range(0, (len(read))):
				#TODO: use named tuples after merge with master branch: read[j].position
				variant_index= vcf_indices.get(read[j].position)
				if variant_index!= None:
					SNPset.append(variant_index)
					score += 1

			#Check for paired end reads and then
			# changed score subtract SNPs covered physically, but not sequenced
			if len(SNPset)!= (SNPset[len(SNPset)-1] - SNPset[0]+1):
				score = score - ((SNPset[len(SNPset)-1] - SNPset[0]+1)-len(SNPset))

			covered_SNPs = tuple(SNPset)

			# TODO: is should be sufficient to only store the index: priorityqueue.push(score,index)
			priorityqueue.push(score,(index,covered_SNPs))
			# TODO: ... also here
			pq_item = tuple([index,covered_SNPs])
			for m in range(0, len(covered_SNPs)):
				SNP_read_map[covered_SNPs[m]].append(pq_item)
		print("Number of skipped reads:", skipped_reads)

		return (priorityqueue,SNP_read_map)



	#TODO Need to assert somewhere that if less Reads than coverage...?
	#TODO insert the bridging .....
	def read_selection(self, pq,SNP_dict, max_coverage):
		'''Selects the best reads out of the readset depending on the assigned score till max_coverage is reached'''

		# coverage over the SNPs
		coverages = CovMonitor(len(self.positions))

		selected_reads = []

		while not pq.is_empty():
			(read,SNPs)=pq._getitem(0)
			pq.pop()

			#Setting coverage of extracted read
			begin = SNPs[0]
			end = SNPs [len(SNPs) - 1] + 1
			if coverages.max_coverage_in_range(begin, end) < max_coverage:
				coverages.add_read(begin, end)
				selected_reads.append(read)

			#Reducing the score of all reads covering the same SNPS as the selected read
			for i in range(0,len(SNPs)):
				SNP_reads= SNP_dict[SNPs[i]]

				j=0
				#if it is the extracted read, then remove it from the list, else decrease the score by 1.
				while j< len(SNP_reads):
					r=SNP_reads[j]

					if r ==(read,SNPs):
						removable_read= (read,SNPs)
						SNP_reads.remove(removable_read)
					else:
						readindex= pq._get_index(r)
						old_score=pq._getscore(readindex)
						newscore= old_score -1
						#TODO need to set a minium for the score ?
						pq.change_score(r,newscore)
						j +=1

		return selected_reads


