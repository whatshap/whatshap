# TODO
# Class Reads where access through SNP position -DONE in SNP MAP
# implement readscore - DONE in combined method to build up SNP MAP which include the readscore
#Heap implemented for storage of Reads - Done

#include the Coverage Monitor...


# Redefine ComponentFinder by using max value and also the rank (not sure)
#implement heuristic
#return Readset in the former representation
#Erase SNPS which are unphasable

import math


from whatshap._core import PyRead as Read
from whatshap._core import PyReadSet as ReadSet
from whatshap._core import PyDPTable as DPTable
from whatshap._core import PyIndexSet as IndexSet
#from whatshap.scripts.whatshap import CoverageMonitor as CovMonitor


#Not possible to import
class CovMonitor:
	'''TODO: This is a most simple, naive implementation. Could do this smarter.'''
	def __init__(self, length):
		self.coverage = [0] * length

	def max_coverage_in_range(self, begin, end):
		return max(self.coverage[begin:end])

	def add_read(self, begin, end):
		for i in range(begin, end):
			self.coverage[i] += 1







class Bestreads:
	#Value is the read
	#parent is the vcf_variant
	def __init__(self, value, parent):
		self.value = value
		self.parent = parent


	#Builds the structure for the variant to Read mapping
	def snp_map_construction(self):
		#set up an empty dictionary with the length of the vcf variants
		d = {}
		for i in range(0, len(self.parent)):
			d[i] = []
		#indexing the reads
		readindices = list(range(len(self.value)))
		#indexing the vcf_variants
		vcf_indices = {position: index for index, position in enumerate(self.parent)}
		skipped_reads = 0
		#if we want to see which SNPs are unphasable need to compute all SNP positions between begin and end position.which may contribute to coverage but not to phasability
		possible_reads = set()

		for index in readindices:
			read = self.value[index]


			if len(read) < 2:
				skipped_reads += 1
				continue

			b_position, b_base, b_allele, b_quality = read[0]
			e_position, e_base, e_allele, e_quality = read[len(read) - 1]
			begin = vcf_indices[b_position]
			#version of whatshap.py by end +1
			#TODO the difference between the 2 methods
			end = vcf_indices[e_position]

			#score for the sorting of best reads
			realscore = 0
			#set containing for each read which variants are covered by the read
			realset = []

			#for the reads and the vcf indices
			for j in range(0, (len(read))):
				for k in range(0, len(vcf_indices)):
					#looking exactly which SNPS are covered in the read
					if (read[j][0]) == (self.parent[k]):
						realset.append(k)
						realscore += 1

			helptupel = [index, realscore, realset]
			for m in range(0, len(realset)):
				d[realset[m]].append(helptupel)

		return d

#Implement this method
#readselect.score_selection(heapstructure,max_coverage)


	def heapcon(self, SNP_map):
		#Storage for the heaps of the SNP variants
		heapstorage=[]

		for i in range(0,len(SNP_map)):
			helpheap=self.create_new_heap(SNP_map[i])
			heapstorage.insert(i,helpheap)
			#actualmax=self.find_max(helpheap)


		#out = [heapstorage,actualmax]
		#TESTING
		#only for testing
		testarray_for_heap_analysis = {}
		testarray_for_heap_analysis[0] = [[1, 2, [0, 1]], [2, 2, [0, 1]], [5, 4, [0, 1, 2, 3]], [6, 3, [0, 1, 2]],[7, 4, [0, 1, 2, 3]]]
		testarray_for_heap_analysis[1] = [[1, 2, [0, 1]], [2, 2, [0, 1]], [3, 2, [1, 2]], [4, 3, [1, 2, 3]], [5, 4, [0, 1, 2, 3]], [6, 3, [0, 1, 2]], [7, 4, [0, 1, 2, 3]]]
		testarray_for_heap_analysis[2] = [[3, 2, [1, 2]], [4, 3, [1, 2, 3]], [5, 4, [0, 1, 2, 3]], [6, 3, [0, 1, 2]], [7, 4, [0, 1, 2, 3]]]
		testarray_for_heap_analysis[3] = [[4, 3, [1, 2, 3]], [5, 4, [0, 1, 2, 3]], [7, 4, [0, 1, 2, 3]]]


		selected_reads=self.score_selection(testarray_for_heap_analysis,5)







		return heapstorage

	#TODO Need to assert somewhere that if less Reads than coverage...?
	def score_selection(self,heap,max_coverage):


		#For testing
		setted_coverage= 4

		# ... and the corresponding coverages along each slice
		coverages = CovMonitor(len(self.parent))

		#Storage of the selected reads
		bestreads=[]

		#current global maximal read
		globalmax = [0,0,[]]
		#current SNP variants which includes the maximum read
		maxindex=[]
		#SNP variant where globalmax is actual the maximum
		maxvalue=-1

		#TODO put it in one method no code copy to method below
		for i in range(0,len(heap)):
			print('i')
			print(i)
			if not self.is_Empty(heap[i]):
				actualmax=self.find_max(heap[i])
				if actualmax[1]> globalmax[1]:
					globalmax=actualmax
					maxindex=globalmax[2]
					maxvalue=i


		begin=globalmax[2][0]
		end= globalmax[2][len(globalmax[2])-1]

		nextiteration=[0]*(len(heap)-1)

		if coverages.max_coverage_in_range(begin,end) < setted_coverage:
			coverages.add_read(begin,end)



			#over all variants which also iniclude this selected read erase the read out of the heap
			for l in maxindex:
				posindex= heap[l].index(globalmax)
				del heap[l][posindex]
				self.max_new_heapify(heap[l],posindex)



			maxindex.remove(maxvalue)
			bestreads.append(globalmax)
			self.Extract_new_Max(heap[maxvalue])
			self.Shift_down(heap[maxvalue])
			nextiteration.insert(maxvalue,heap[maxvalue])
			heap[maxvalue]=[]


		#1 iteration already  done
		numbervariants=len(heap)-1
		#search exactly number of variants time for the best read so that each SNP should be  covered once.
		while numbervariants!= 0:
			(heap,nextiteration,bestreads,coverages)=self.get_max_read(setted_coverage,heap,nextiteration,bestreads,coverages)
			#print('while')
			#print(heap)
			#print(nextiteration)
			#print(bestreads)
			#print(coverages.max_coverage_in_range(0,len(self.parent)))
			numbervariants -=1


		return bestreads


	def get_max_read(self,max_coverage,oldheap,newheap,bestreads,cmonitor):
		#current global maximal read
		globalmax = [0,0,[]]
		#current SNP variants which includes the maximum read
		maxindex=[]
		#SNP variant where globalmax is actual the maximum
		maxvalue=-1

		for i in range(0,len(oldheap)):
			print('i')
			print(i)
			if not self.is_Empty(oldheap[i]):
				actualmax=self.find_max(oldheap[i])
				if actualmax[1]> globalmax[1]:
					#print('max found')
					globalmax=actualmax
					maxindex=globalmax[2]
					maxvalue=i
			else:
				continue

		if globalmax!=[0,0,[]]:

			begin=globalmax[2][0]
			end= globalmax[2][len(globalmax[2])-1]

			if cmonitor.max_coverage_in_range(begin,end) < max_coverage:
				cmonitor.add_read(begin,end)

				for l in maxindex:

					if not self.is_Empty(oldheap[l]):
						posindex= oldheap[l].index(globalmax)
						del oldheap[l][posindex]
						self.max_new_heapify(oldheap[l],posindex)

				maxindex.remove(maxvalue)
				bestreads.append(globalmax)
				self.Extract_new_Max(oldheap[maxvalue])
				self.Shift_down(oldheap[maxvalue])
				newheap.insert(maxvalue,oldheap[maxvalue])
				oldheap[maxvalue]=[]

		out=[oldheap,newheap,bestreads,cmonitor]
		return out






#math.ceil  aufrunden
#math.floor abrunden

	#maybe also coded by bitshifting
	def hparent(self,i,h):
		return (math.ceil(i/2)-1)
	#maybe also coded by bitshifting
	def hleft(self,i,h):
		return (math.ceil(2*i)+1)
	#maybe also coded by bitshifting
	def hright(self,i,h):
		return (math.ceil(2*i)+2)


	def extract_Element(self,H,read):
		H.remove(read)
		return 0




	def create_new_heap(self,A):
		heapsize=self.size(A)
		index = math.floor((heapsize/2))-1
		while index >= 0:
			self.max_new_heapify(A,index)
			index -=1
		return A





	#create a heap out of given array of elements
	def max_new_heapify(self, A,i):
		l=self.hleft(i,A)
		r=self.hright(i,A)
		#TODO STILL NOT SURE IF INDEX IS REALLY RIGHT....
		if l<= (self.size(A)-1) and A[l][1]>A[i][1]:
			largest=l
		else:
			largest=i
			#NOT sure if small or like (in book smaller equal)
		if r < self.size(A) and A[r][1] >A[largest][1]:
			largest = r
		if largest != i:
			help= A[i]
			A[i]=A[largest]
			A[largest]=help
			self.max_new_heapify(A,largest)




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
	def size(self,h):
		size = len(h)
		return size

	#returns true if the heap is empty, false otherwise.
	def is_Empty(self,H):
		if len(H)==0:
			return True
		else:
			return False

	#TODO Throw Error
	def Extract_new_Max(self,H):
		if self.size(H)<1:
			print('Error heap underflow')
		max= self.find_max(H)
		H[0]=H[self.size(H)-1]
		del H[len(H)-1]
		self.max_new_heapify(H,0)

		return max




	#decreases the scores of all reads in the heap
	def Shift_down(self,H):
		#dont want negative scores
		for i in range(0,len(H)):
			if H[i][1]!= 0 :
				H[i][1]-=1
			else :
				print('Score already below 0')
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