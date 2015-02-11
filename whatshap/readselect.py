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
		la = PriorityQueue()
		print('LALALALALa')
		print(la)
		self.priorityqueue_construction(la)

	def priorityqueue_construction(self, priorityqueue):
		'''Constructiong of the priority queue for the readset, so that each read is in the priority queue and sorted by their score'''
		d = {}
		for i in range(0, len(self.positions)):
			d[i] = []

		#indexing the reads
		readindices = list(range(len(self.readset)))
		#dictionary to map the SNP positions to an index
		vcf_indices = {position: index for index, position in enumerate(self.positions)}

		skipped_reads = 0
		#if we want to see which SNPs are unphasable need to compute all SNP positions between begin and end position.which may contribute to coverage but not to phasability

		for index in readindices:
			read = self.readset[index]

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
				#looking exactly which variants are covered in the read
				variant_index= vcf_indices.get(read[j][0])
				if variant_index!= None:
					SNPset.append(variant_index)
					score += 1
			#Check for paired end reads and then
			# changed score subtract SNPs covered physically, but not sequenced
			if len(SNPset)!= (SNPset[len(SNPset)-1] - SNPset[0]+1):
				score = score - ((SNPset[len(SNPset)-1] - SNPset[0]+1)-len(SNPset))

			covered_SNPs = tuple(SNPset)
			#print('SNP list')
			#print(SNP_list)
			priorityqueue.push(score,(index,covered_SNPs))
			helptupel = tuple([index, score, covered_SNPs])
			for m in range(0, len(covered_SNPs)):
				d[covered_SNPs[m]].append(helptupel)



#Only need to store which SNPs ar covered by which read



	def snp_map_construction(self):
		'''Builds the structure for the variant to Read mapping.
        TODO: say what "the structure" is exactly
        '''

		#set up an empty dictionary with the length of the vcf variants
		d = {}
		for i in range(0, len(self.positions)):
			d[i] = []


		#indexing the reads
		readindices = list(range(len(self.readset)))
		#indexing the vcf_variants
		vcf_indices = {position: index for index, position in enumerate(self.positions)}
		#print('vcf indicies')
		#print(vcf_indices)
		skipped_reads = 0
		#if we want to see which SNPs are unphasable need to compute all SNP positions between begin and end position.which may contribute to coverage but not to phasability

		for index in readindices:
			read = self.readset[index]

			if len(read) < 2:
				skipped_reads += 1
				continue

			#score for the sorting of best reads
			realscore = 0
			#set containing for each read which variants are covered by the read
			realset = []

			#for the reads and the vcf indices
			for j in range(0, (len(read))):
				#TODO: don't iterate, query dictionary vcf_indices[position]
				for k in range(0, len(vcf_indices)):
					#TODO: use named tuples after merge with master branch: read[j].position
					#looking exactly which SNPS are covered in the read
					if (read[j][0]) == (self.positions[k]):
						realset.append(k)
						realscore += 1
			# TODO: subtract SNPs covered physically, but not sequenced

			helptupel = [index, realscore, realset]
			for m in range(0, len(realset)):
				d[realset[m]].append(helptupel)

		return d

	#Size of the heapstorage or mor precise the number of variants which were still active...
	def heapstoragesize(self, heapstorage):
		size = 0
		for i in range(0, len(heapstorage)):
			if not self.is_Empty(heapstorage[i]):
				size += 1
		return size


	def heapcon(self, SNP_map, max_coverage):
		#Storage for the heaps of the SNP variants
		heapstorage = []
		unphasable_SNP_positions = 0

		for i in range(0, len(SNP_map)):


			helpheap = self.create_new_heap(SNP_map[i])
			heapstorage.insert(i, helpheap)
			#actualmax=self.find_max(helpheap)
			if self.is_Empty(SNP_map[i]):
				unphasable_SNP_positions += 1






		#for i in range(0,len(testarray_for_heap_analysis)):
		#	helpheap=self.create_new_heap(testarray_for_heap_analysis[i])
		#	heapstorage.insert(i,helpheap)

		selected_reads2 = self.score_selection(heapstorage, max_coverage)

		print('Found %d  unphasable SNPs ' % unphasable_SNP_positions)
		return selected_reads2

	#TODO Need to assert somewhere that if less Reads than coverage...?
	def score_selection(self, heap, max_coverage):

		# ... and the corresponding coverages along each slice
		coverages = CovMonitor(len(self.positions))
		#print('COVERAGE MONITOR INITIALIZED')


		#Storage of the selected reads
		bestreads = []
		nextiteration = [[] for i in range(len(heap))]
		bestreads2 = []

		Going_on = True

		while Going_on:

			while self.heapstoragesize(heap) != 0 and Going_on:
				(heap, nextiteration, bestreads, coverages, Going_on) = self.get_max_read(max_coverage, heap,
																						  nextiteration, bestreads,
																						  coverages, Going_on)

			#TODO between here the union find algorithm
			#Maybe distinguish between real covered SNPs and SNPS that are only covered but are actually unphasable

			while self.heapstoragesize(nextiteration) != 0 and Going_on:
				(nextiteration, heap, bestreads2, coverages, Going_on) = self.get_max_read(max_coverage, nextiteration,
																						   heap, bestreads2, coverages,
																						   Going_on)

		selected = bestreads + bestreads2

		return selected


	def get_max_read(self, max_coverage, oldheap, newheap, bestreads, cmonitor, Going_on):


		#current global maximal read
		globalmax = [0, -1, []]
		#current SNP variants which includes the maximum read
		maxindex = []
		#SNP variant where globalmax is actual the maximum
		maxvalue = -1

		for i in range(0, len(oldheap)):
			if not self.is_Empty(oldheap[i]):
				actualmax = self.find_max(oldheap[i])


				#Set globalmax score to -1 because else the score 0 is not possible to see.
				if actualmax[1] > globalmax[1]:
					#print('max found')
					globalmax = actualmax
					maxindex = globalmax[2]

					maxvalue = i

			else:
				continue

		#TODO if globalmax is found otherwise ?????
		if globalmax != [0, -1, []]:

			begin = globalmax[2][0]
			end = globalmax[2][len(globalmax[2]) - 1] + 1

			if cmonitor.max_coverage_in_range(begin, end) < max_coverage:
				cmonitor.add_read(begin, end)
				#print('read added')

				bestreads.append(globalmax)

				#Go over all variants which are covered by this read
				#not -1 because this would accesss the last element of the list
				posindex = -10
				secposindex = -10
				for l in maxindex:
					if not self.is_Empty(oldheap[l]):
						#Go over the heaps of each variant
						for j in range(0, len(oldheap[l])):
							#Compare if they are the same Read by the Read indices
							if oldheap[l][j][0] == globalmax[0]:
								posindex = j
							#print(posindex)
						#only delete the read if found (by paired end reads not essentially in the read also in the heap of the covered SNP variant)
						if posindex != -1:
							del oldheap[l][posindex]

						self.Shift_down(oldheap[l])
						self.max_new_heapify(oldheap[l], posindex)

					#Same for newheap because covered SNP could already be moved in the second heapstructure

					else:
						#if l <= len(newheap):
						#print('ZU klein ')
						if not self.is_Empty(newheap[l]):
							for j in range(0, len(newheap[l])):
								if newheap[l][j][0] == globalmax[0]:
									secposindex = j
								#print(secposindex)
							if secposindex != -1:
								del newheap[l][secposindex]
							self.Shift_down(newheap[l])
							self.max_new_heapify(newheap[l], secposindex)

				maxindex.remove(maxvalue)
				#only need to to if in the heap something is left... after removing the globalmax
				if not self.is_Empty(oldheap[maxvalue]):
					#TODO is not needed the same is already done by posindex maybe use this method..
					#self.Extract_new_Max(oldheap[maxvalue])
					newheap.pop(maxvalue)
					newheap.insert(maxvalue, oldheap[maxvalue])
					oldheap[maxvalue] = []

		if cmonitor.max_coverage_in_range(0, len(oldheap)) == max_coverage:
			Going_on = False

		out = [oldheap, newheap, bestreads, cmonitor, Going_on]

		return out


	#math.ceil  aufrunden
	#math.floor abrunden

	#maybe also coded by bitshifting
	def hparent(self, i, h):
		return (math.ceil(i / 2) - 1)

	#maybe also coded by bitshifting
	def hleft(self, i, h):
		return (math.ceil(2 * i) + 1)

	#maybe also coded by bitshifting
	def hright(self, i, h):
		return (math.ceil(2 * i) + 2)


	def extract_Element(self, H, read):
		H.remove(read)
		return 0


	def create_new_heap(self, A):
		heapsize = self.size(A)
		index = math.floor((heapsize / 2)) - 1
		while index >= 0:
			self.max_new_heapify(A, index)
			index -= 1
		return A


	#create a heap out of given array of elements
	def max_new_heapify(self, A, i):
		l = self.hleft(i, A)
		r = self.hright(i, A)
		#TODO STILL NOT SURE IF INDEX IS REALLY RIGHT....
		if l <= (self.size(A) - 1) and A[l][1] > A[i][1]:
			largest = l
		else:
			largest = i
		#NOT sure if small or like (in book smaller equal)
		if r < self.size(A) and A[r][1] > A[largest][1]:
			largest = r
		if largest != i:
			help = A[i]
			A[i] = A[largest]
			A[largest] = help
			self.max_new_heapify(A, largest)


	#TODO include here if heap is empty needed ?
	#find the maximum item of a max-heap
	def find_max(self, H):
		h = H[0]
		return h

	'''
    #removing the root node of a max-heap, AND decrease the score of the other reads included in the SNPS covered by this read by 1
    def delete_max(self, H):
        max=self.Extract_new_Max(H)
        #TODO MAYBE INCLUDED INTO EXTRACT MAX
        #covered_variants = self.Variants_of_Max
        #for k in covered_variant.getReads :
        #	deccrease_key (k)
        #	heapify Reads (index ....)
        h = 0
        return h
    '''

	#return the number of items in the heap.
	def size(self, h):
		size = len(h)
		return size

	#returns true if the heap is empty, false otherwise.
	def is_Empty(self, H):
		if len(H) == 0:
			return True
		else:
			return False

	#TODO Throw Error
	def Extract_new_Max(self, H):
		if self.size(H) < 1:
			print('Error heap underflow')
		max = self.find_max(H)
		H[0] = H[self.size(H) - 1]
		del H[len(H) - 1]
		self.max_new_heapify(H, 0)

		return max


	#decreases the scores of all reads in the heap
	def Shift_down(self, H):
		#dont want negative scores
		for i in range(0, len(H)):
			if H[i][1] != 0:
				H[i][1] -= 1
			#else :
			#print('Score already below 0')
			#Nothing to be done only for information...
			#TODO better solution .... Remove node out of HEAP ?

		return 0


##############################################

#code from graph.py to modify Componentfinder for the read selection purpose

class Node:
	def __init__(self, value, parent):
		self.value = value
		self.parent = parent

	def __repr__(self):
		return "Node(value={}, parent={})".format(self.value, self.parent)


class ComponentFinder2:
	"""
    Find connected components. A ComponentFinder is initialized with a list of
    values. These are initially partitioned such that each value is in a
    separate set. By calling merge(x, y), the two sets containing values x and
    y are merged. Calling find(x) returns a "representative" value of the set
    that value x is in. x and y are in the same set iff find(x) == find(y).
    The representative is always the minimum value of the set.

    # TODO:need to have the maximum value of the set because we search the max of the readscore.



    This implements a variant of the Union-Find algorithm, but without the
    "union by rank" strategy since we want the smallest node to be the
    representative. It could perhaps be optimized, but this function is not
    the current bottleneck.
    """

	def __init__(self, values):
		self.nodes = {x: Node(x, None) for x in values}

	def merge(self, x, y):
		assert x != y
		x_root = self._find_node(x)
		y_root = self._find_node(y)

		if x_root is y_root:
			return

		# Merge while making sure that the node with the smaller value is the
		# new parent.
		if x_root.value < y_root.value:
			y_root.parent = x_root
		else:
			x_root.parent = y_root

	def _find_node(self, x):
		node = root = self.nodes[x]
		while root.parent is not None:
			root = root.parent

		# compression path
		while node.parent is not None:
			node.parent, node = root, node.parent
		return root

	def find(self, x):
		"""
        Return which component x belongs to, identified by the smallest value.
        """
		return self._find_node(x).value

	def print(self):
		for x in sorted(self.nodes):
			print(x, ':', self.nodes[x], 'is represented by', self._find_node(x))