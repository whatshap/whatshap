import math


class PriorityQueue:

	def __init__(self):
		'''initializes a priority queue with an empty list and an empty dictionary'''
		self.positions={}
		self.heap=[]


	def push(self, score, item):
		'''Add item with given score to the heap.'''

		pass

	def pop(self):
		'''Removes the item with largest score and returns it as (score, item).'''
		pass

	def change_score(self, item, new_score):
		'''Changes the score of the given item.'''
		#TODO: to implement that, one could use a dictionary mapping items to its position in the heap
		pass

	def __len__(self):
		return 0


	#Indexes of parent, left and right child maybe also coded by bitshifting

	def hparent(self,i,h):
		''''returns the index of the parent node in the heap'''
		return (math.ceil(i/2)-1)

	def hleft(self,i,h):
		'''return the index of the left child in the heap'''
		return (math.ceil(2*i)+1)

	def hright(self,i,h):
		'''return the index of the right child in the heap'''
		return (math.ceil(2*i)+2)