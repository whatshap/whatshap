
class PriorityQueue:
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
