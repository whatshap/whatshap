import math


class PriorityQueue:
	def __init__(self):
		'''initializes a priority queue with an empty list and an empty dictionary'''
		self.positions = {}
		self.heap = []


	def push(self, score, item):
		'''Add item with given score to the heap.'''
		# Item stored with score in a tupel where the score is always the first position..
		# insert item in the last position of the heap and then change the position till it suits
		newindex = len(self.heap) + 1
		self.heap.insert(newindex, (score, item))
		self.positions[item] = newindex
		self.parent_check_and_swap(newindex)

		#if self.parent_check_and_swap(newindex):
		#	continue
		#else:
		#	print('Pushing does not work')


		#TODO not sure if heapify in this case needed
		#self.heapify(newindex)



	def swap(self, f_index, s_index):
		'''swaps the position of the nodes in the priority queue from first index(f_index)and second index(s_index)'''
		firstitem = self.heap[f_index]
		firstpos = self.positions[self.getitem(f_index)]
		self.heap[f_index] = self.heap[s_index]
		self.positions[self.getitem(f_index)] = self.positions[self.getitem(s_index)]
		self.heap[s_index] = firstitem
		self.positions[self.getitem(s_index)] = firstpos


	def parent_check_and_swap(self, index):
		''' Check if score of item at given index is higher than the score of the parent and swaps the nodes till priority
		queue property is restored'''
		parentindex = self.hparent(index)
		if self.getscore(parentindex) < self.getscore(index):
			self.swap(parentindex, index)
			# Look if that works with the Index after swapping
			self.parent_check_and_swap(parentindex)
		#else:
		#	return

	def child_check_and_swap(self,index):



	def pop(self):
		'''Removes the item with largest score and returns it as (score, item).'''
		largest_element = self.extract_max()
		largest_element = self.heap[0]
		# Not sure if this works (TODO TEST)
		self.positions.pop(largest_element)

		return self.heap[0]
		pass


	def change_score(self, item, new_score):
		'''Changes the score of the given item.'''
		position=self.positions[item]
		old_score= self.getscore(position)
		self.heap.insert(position,(new_score,item))
		if old_score<new_score:
			self.parent_check_and_swap(position)
		else:
			self.child_check_and_swap(position)

				# Differentiate between higher or lower values in order to compare to the parent or to the children....

		# TODO: to implement that, one could use a dictionary mapping items to its position in the heap



# TODO do not know if this method is required.
	def getscore(self, index):
		'''Return the score of the element at this index in the heap'''
		(score, item) = self.heap[index]
		return score

	def getitem(self, index):
		'''Return the score of the element at this index in the heap'''
		(score, item) = self.heap[index]
		return item


	def __len__(self):
		'''Length of  Priority Queue is equal to the length of the stored heap'''
		return len(self.heap)


	def isEmpty(self):
		'''Return if sctual Priority Queue is Empty'''
		if len(self) == 0:
			return True
		else:
			return False


	def heapify(self, index):
		'''Restores the max priority queue property starting at the given index'''
		l = self.hleft(index)
		r = self.hright(index)

		if l <= (len(self) - 1) and self.getscore(l) > self.getscore(index):
			largest = l
		else:
			largest = index
		# NOT sure if small or like (in book smaller equal)
		if r < len(self) and self.getscore(r) > self.getscore(largest):
			largest = r
		if largest != index:
			#change the index and largest node in the heap
			help = self.heap[index]
			helppos = self.positions[index]
			self.heap[index] = self.heap[largest]
			self.positions[index] = self.positions[largest]
			self.heap[largest] = help
			self.positions[largest] = helppos

			self.heapify(largest)


	def extract_max(self):
		'''Gives maximal item in the priority queue , erases it from the queue and updates the scores and indices'''


	def compare(self, item1, item2):


	# Indexes of parent, left and right child maybe also coded by bitshifting

	def hparent(self, index):
		''''returns the index of the parent node in the heap depending on the given index'''
		return (math.ceil(index / 2) - 1)


	def hleft(self, index):
		'''return the index of the left child in the heap depending on the given index'''
		return (math.ceil(2 * index) + 1)


	def hright(self, index):
		'''return the index of the right child in the heap depending on the given index'''
		return (math.ceil(2 * index) + 2)


	#def extract_Element(self,H,read):
	#	H.remove(read)
	#	return 0



	#Not needed if we push each element for itself...
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





