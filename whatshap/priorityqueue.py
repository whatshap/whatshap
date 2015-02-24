import math

#Important to use tuples instead of lists of SNPs because llist are unhasable for the dictionary

# TODO: prefix with underscore: getscore, getitem, get_index, hparent, hleft, hright, swaps

# TODO: isEmpty --> is_empty

# TODO: in general, see https://www.python.org/dev/peps/pep-0008

class PriorityQueue:
	def __init__(self):
		'''initializes a priority queue with an empty list and an empty dictionary'''
		self.positions = {}
		self.heap = []

	def push(self, score, item):
		'''Add item with given score to the heap.'''
		# Item stored with score in a tupel where the score is always the first position..
		newindex = len(self.heap)
		self.heap.insert(newindex, (score, item))
		self.positions[item] = newindex
		self._sift_up(newindex)



	def _swap(self, f_index, s_index):
		'''swaps the position of the nodes in the priority queue from first index(f_index)and second index(s_index)'''

		firstitem = self.heap[f_index]
		firstpos = self.positions[self._get_item(f_index)]

		seconditem = self.heap[s_index]
		secpos= self.positions[self._get_item(s_index)]

		self.positions[self._get_item(f_index)] =secpos
		self.positions[self._get_item(s_index)] =firstpos

		self.heap[f_index] = seconditem
		self.heap[s_index] = firstitem


	def _sift_up(self, index):
		''' recursive method to check if score of item at given index is higher than the score of the parent
		and swaps the nodes till priorityqueue property is restored'''
		parentindex = self._hparent(index)
		assert parentindex != index
		if parentindex >=0 :
			if self._get_score(parentindex) < self._get_score(index) :
				self._swap(parentindex, index)
				self._sift_up(parentindex)


	def _sift_down(self,index):
		'''Check if score of item at given index is lower than the score of its children,
		 therefore need to swap position with its children'''
		rchildindex=self._hright(index)
		lchildindex=self._hleft(index)
		assert rchildindex != index
		assert lchildindex != index
		#if both children are in the heap
		#only need to know if right child exists then automatically also the left child exists
		if rchildindex<len(self.heap):
			rchildscore= self._get_score(rchildindex)
			lchildscore= self._get_score(lchildindex)
			#returns true if right child has an higher score as left child
			score_comparison= rchildscore >lchildscore

			if score_comparison and rchildindex<len(self.heap):
				'''right child has higher score than the actual index '''
				self._swap(rchildindex,index)
				self._sift_down(rchildindex)
			else:
				if lchildscore>self._get_score(index) and lchildindex<len(self.heap):
					'''left child is smaller than right child but still bigger than index'''
					self._swap(lchildindex,index)
					self._sift_down(lchildindex)
		else:
			if lchildindex<len(self.heap) and self._get_score(lchildindex)>self._get_score(index):
				self._swap(lchildindex,index)
				self._sift_down(lchildindex)




	def pop(self):
		'''Removes the item with largest score and returns it as tupel of (score, item).'''

		#looks if heap is only one element then no need restore heap property
		if len(self.heap)>1:
			#remember last element, max element and their  items
			last_element=self.heap[len(self.heap)-1]

			item_latest=self._get_item(len(self.heap)-1)
			item_max=self._get_item(0)

			max_element =self.heap.pop(0)

			self.heap.pop(len(self.heap)-1)
			self.heap.insert(0,last_element)

			self.positions[item_latest]= 0
			self.positions.pop(item_max)

			self._sift_down(0)

		else:
			self.positions.pop(self._get_item(0))
			max_element=self.heap.pop(0)

		return max_element




	def change_score(self, item, new_score):
		'''Changes the score of the given item to the new assigned score.'''
		position=self.positions[item]
		old_score= self.get_score_by_item(item)
		self.heap.pop(position)

		self.heap.insert(position,(new_score,item))

		# Differentiate between increasing and decreasing score
		if old_score<new_score:
			self._sift_up(position)
		else:
			self._sift_down(position)	


	def _get_score(self, index):
		'''Return the score of the element at this index in the heap'''
		(score, item) = self.heap[index]
		return score

	def _get_item(self, index):
		'''Return the score of the element at this index in the heap'''
		(score, item) = self.heap[index]
		return item

	def _get_index (self,item):
		'''returns intern index of the searched item'''
		return self.positions[item]

	def get_score_by_item(self,item):
		'''returns actual score of the given item '''
		(score,item) = self.heap[self._get_index(item)]
		return score



	def __len__(self):
		'''Length of  Priority Queue is equal to the length of the stored heap'''
		return len(self.heap)


	def is_empty(self):
		'''Return if actual Priority Queue is Empty'''
		if len(self) == 0:
			return True
		else:
			return False

	#TODO maybe code it via bitshifting

	def _hparent(self, index):
		''''returns the index of the parent node in the heap depending on the given index'''
		return (math.ceil(index / 2) - 1)


	def _hleft(self, index):
		'''return the index of the left child in the heap depending on the given index'''
		return (math.ceil(2 * index) + 1)


	def _hright(self, index):
		'''return the index of the right child in the heap depending on the given index'''
		return (math.ceil(2 * index) + 2)


