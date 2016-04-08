from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.unordered_map cimport unordered_map
from libcpp cimport bool

import math


cdef bool _vector_score_lower(priority_type_ptr first, priority_type_ptr second):
	for i in range(min(first[0].size(), second[0].size())):
		if first[0][i] < second[0][i]:
			return True
		if first[0][i] > second[0][i]:
			return False
	if first[0].size() < second[0].size(): 
		return True
	else:
		return False


cdef priority_type_ptr _pyscore_to_vector(score):
	cdef priority_type_ptr result = new priority_type()
	if isinstance(score,int):
		result[0].push_back(score)
	else:
		try:
			for i in score:
				if not isinstance(i,int):
					raise
				result[0].push_back(i)
		except:
			raise ValueError('Score parameter must be either int, or an iterable object yielding ints')
	return result


cdef int _parent(int index):
	''''returns the index of the parent node in the heap depending on the given index'''
	return (index - 1) // 2


cdef int _left_child(int index):
	'''return the index of the left child in the heap depending on the given index'''
	return (2 * index) + 1


cdef int _right_child(int index):
	'''return the index of the right child in the heap depending on the given index'''
	return (2 * index) + 2


# TODO: in general, see https://www.python.org/dev/peps/pep-0008
cdef class PriorityQueue:
	def __dealloc__(self):
		for i in range(self.heap.size()):
			del self.heap[i].first
	
	def push(self, score, int item):
		'''Add item with given score to the heap.'''
		# Item stored with score in a tupel where the score is always the first position..
		cdef priority_type_ptr c_score = _pyscore_to_vector(score)
		self.c_push(c_score, item)

	cdef void c_push(self, priority_type_ptr score, int item):
		'''Pointer ownership is transferred to priorityqueue.'''
		newindex = self.heap.size()
		cdef queue_entry_type entry
		entry.first = score
		entry.second = item 
		self.heap.push_back(entry)
		self.positions[item] = newindex
		self._sift_up(newindex)

	cdef void _swap(self, int index1, int index2):
		'''swaps the position of the nodes in the priority queue from first index(index1)and second index(index2)'''

		cdef queue_entry_type entry1 = self.heap[index1]
		cdef int pos1 = self.positions[entry1.second]

		cdef queue_entry_type entry2 = self.heap[index2]
		cdef int pos2 = self.positions[entry2.second]

		self.positions[entry1.second] = pos2
		self.positions[entry2.second] = pos1

		self.heap[index1] = entry2
		self.heap[index2] = entry1

	cdef bool _score_lower(self, int index1, int index2):
		cdef queue_entry_type entry1 = self.heap[index1]
		cdef queue_entry_type entry2 = self.heap[index2]
		return _vector_score_lower(entry1.first, entry2.first)

	cdef void _sift_up(self, int index):
		''' recursive method to check if score of item at given index is higher than the score of the parent
		and swaps the nodes till priorityqueue property is restored'''
		#print(str(self), '_sift_up', index)
		cdef int parentindex = _parent(index)
		assert parentindex != index
		if parentindex >= 0:
			if self._score_lower(parentindex, index):
				self._swap(parentindex, index)
				self._sift_up(parentindex)


	cdef void _sift_down(self, int index):
		'''Check if score of item at given index is lower than the score of its children,
		 therefore need to swap position with its children'''
		cdef int rchildindex = _right_child(index)
		cdef int lchildindex = _left_child(index)
		assert rchildindex != index
		assert lchildindex != index
		#if both children are in the heap
		#only need to know if right child exists then automatically also the left child exists
		if rchildindex < self.heap.size():
			if self._score_lower(lchildindex, rchildindex):
				if self._score_lower(index, rchildindex):
					self._swap(rchildindex, index)
					self._sift_down(rchildindex)
			else:
				if self._score_lower(index, lchildindex):
					self._swap(lchildindex, index)
					self._sift_down(lchildindex)
		elif lchildindex < self.heap.size():
			if self._score_lower(index, lchildindex):
				self._swap(lchildindex, index)
				self._sift_down(lchildindex)

	def pop(self):
		'''Removes the item with largest score and returns it as tupel of (score, item).'''
		cdef queue_entry_type entry = self.c_pop()
		cdef priority_type_ptr score = entry.first
		cdef item_type item = entry.second
		if score[0].size() == 1:
			result = (score[0][0], item)
		else:
			result = tuple(score[0]), item
		del score
		return result

	cdef queue_entry_type c_pop(self):
		'''Removes the item with largest score and returns it as tupel of (score, item).'''
		if self.heap.size() == 0:
			raise IndexError('PriorityQueue empty.')
		cdef queue_entry_type last_entry = self.heap[self.heap.size()-1]
		cdef queue_entry_type first_entry = self.heap[0]
		#looks if heap is only one element then no need restore heap property
		if self.heap.size() == 1:
			self.positions.erase(first_entry.second)
			self.heap.pop_back()
		else:
			self.heap[0] = last_entry
			self.heap.pop_back()

			self.positions[last_entry.second]= 0
			self.positions.erase(first_entry.second)

			self._sift_down(0)
		return first_entry

	def change_score(self, item_type item, new_score):
		'''Changes the score of the given item to the new assigned score.'''
		cdef priority_type_ptr c_new_score = _pyscore_to_vector(new_score)
		self.c_change_score(item, c_new_score)

	cdef void c_change_score(self, item_type item, priority_type_ptr c_new_score):
		'''Changes the score of the given item to the new assigned score.
		Ownership of c_new_score goes to priorityqueue.'''
		position = self.positions[item]

		c_old_score = self.heap[position].first
		self.heap[position].first = c_new_score

		# Differentiate between increasing and decreasing score
		if _vector_score_lower(c_old_score, c_new_score):
			self._sift_up(position)
		else:
			self._sift_down(position)

		del c_old_score

	def get_score_by_item(self, item):
		'''returns actual score of the given item or None if item is not in the heap '''
		cdef priority_type_ptr score = self.c_get_score_by_item(item)
		if score == NULL:
			return None
		elif score[0].size() == 1:
			return score[0][0]
		else:
			return tuple(score[0])

	cdef priority_type_ptr c_get_score_by_item(self, item_type item):
		'''Returns score of the given item or NULL if item is not in the heap. 
		Pointer ownership stays with priorityqueue.'''
		cdef unordered_map[item_type,int].iterator it = self.positions.find(item)
		if it == self.positions.end():
			return NULL
		cdef queue_entry_type entry = self.heap[self.positions[item]]
		return entry.first
		cdef priority_type_ptr score = entry.first
		if score[0].size() == 1:
			return score[0][0]
		else:
			return tuple(score[0])

	def __len__(self):
		'''Length of  Priority Queue is equal to the length of the stored heap'''
		return self.heap.size()

	cdef int size(self):
		return self.heap.size()

	def is_empty(self):
		'''Return if actual Priority Queue is Empty'''
		return self.heap.size() == 0

	cdef bool c_is_empty(self):
		return self.heap.size() == 0
