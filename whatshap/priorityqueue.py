import math

#TODO : Important to use tuples instead of lists of SNPs because llist are unhasable for the dictionary

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
		self.parent_check_and_swap(newindex)
		#print('Push item')
		#print(item)
		if item[0]==62:
			print('&§ found')
			print(self.positions)
		#TODO Need to look to children not necessary if only pushed to parent.
		#self.child_check_and_swap(newindex)

	def get_index (self,item):
		#print('Positions')
		#print(self.positions)
		#print('Length of Positions')
		#print(len(self.positions))
		return self.positions[item]


	def swap(self, f_index, s_index):
		'''swaps the position of the nodes in the priority queue from first index(f_index)and second index(s_index)'''

		firstitem = self.heap[f_index]
		firstpos = self.positions[self.getitem(f_index)]

		seconditem = self.heap[s_index]
		secpos= self.positions[self.getitem(s_index)]

		self.positions[self.getitem(f_index)] =secpos
		self.positions[self.getitem(s_index)] =firstpos

		self.heap[f_index] = seconditem
		self.heap[s_index] = firstitem


	def parent_check_and_swap(self, index):
		''' Check if score of item at given index is higher than the score of the parent and swaps the nodes till priority
		queue property is restored'''
		parentindex = self.hparent(index)
		if parentindex!= index and parentindex >=0 :
			if self.getscore(parentindex) < self.getscore(index) :

				self.swap(parentindex, index)
				self.parent_check_and_swap(parentindex)


	def child_check_and_swap(self,index):
		'''Check if score of item at given index is lower than the score of its children, therefore need to swap position with
		children'''
		rchildindex=self.hright(index)
		lchildindex=self.hleft(index)
		if rchildindex!= index and rchildindex < len(self.heap):
			if self.getscore(rchildindex)> self.getscore(index):
				self.swap(rchildindex,index)
				self.child_check_and_swap(rchildindex)
		if lchildindex != index and lchildindex < len(self.heap):
			if self.getscore(lchildindex)> self.getscore(index):
				self.swap(lchildindex,index)
				self.child_check_and_swap(lchildindex)



	def pop(self):
		'''Removes the item with largest score and returns it as (score, item).'''
		#Gives maximal item in the priority queue , erases it from the queue and updates the scores and indices'''

		#looks if heap is only one element then no need restore heap property
		if len(self.heap)>1:
			#remember last element, max element and their  items
			last_element=self.heap[len(self.heap)-1]

			item_latest=self.getitem(len(self.heap)-1)
			item_max=self.getitem(0)

			max_element =self.heap.pop(0)

			self.heap.pop(len(self.heap)-1)
			self.heap.insert(0,last_element)

			self.positions[item_latest]= 0
			self.positions.pop(item_max)

			self.child_check_and_swap(0)

		else:
			self.positions.pop(self.getitem(0))
			max_element=self.heap.pop(0)

		return max_element




	def change_score(self, item, new_score):
		'''Changes the score of the given item.'''
		position=self.positions[item]
		old_score= self.getscore(position)
		self.heap.pop(position)
		self.heap.insert(position,(new_score,item))

		# Differentiate between increasing and decreasing score
		if old_score<new_score:
			self.parent_check_and_swap(position)
		else:
			self.child_check_and_swap(position)


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
		'''Return if actual Priority Queue is Empty'''
		if len(self) == 0:
			return True
		else:
			return False

	#TODO maybe code it via bitshifting
	def hparent(self, index):
		''''returns the index of the parent node in the heap depending on the given index'''
		return (math.ceil(index / 2) - 1)


	def hleft(self, index):
		'''return the index of the left child in the heap depending on the given index'''
		return (math.ceil(2 * index) + 1)


	def hright(self, index):
		'''return the index of the right child in the heap depending on the given index'''
		return (math.ceil(2 * index) + 2)


