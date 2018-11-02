from collections import defaultdict
import sys
from .graph import ComponentFinder
from .core import Read, ReadSet
import logging
import pysam
import itertools

logger = logging.getLogger(__name__)

class ErrorsPerColumn:
	"""
	Given a ReadSet which represents the transformed allele matrix, determine the number of corrections 
	in each column needed to make it consistent with the true parition (encoded in the readnames).
	"""
	def __init__(self, reads, ploidy):
		"""
		reads -- ReadSet with reads to be considered
		position_to_component -- dict mapping variant positions to connected components
		ploidy -- ploidy of the data
		"""
		# currently considered positions
		self._position = None
		# ploidy of the data
		self._ploidy = ploidy
		# currently considered window
		self._window = 0
		# set of reads covering current position
		self._current_column = None
		# current positions
		self._positions = []
		# column costs
		self._column_costs = []
		# normalized column costs
		self._normalized_column_costs = []

		# iterate over all positions in the readset
		all_positions = reads.get_positions()
		for i,position in enumerate(all_positions):
			self._current_column = []
			# get current position
			self._position = all_positions[i]
			# get all reads that cover the current position
			for j,read in enumerate(reads):
				add_read = True
				if self._position in read:
					self._current_column.append((j,read))
			if len(self._current_column) == 0:
				continue
			# TODO count errors per column
			self._compute_costs()
			self._window += 1


	def get_column_costs(self):
		"""
		Return the column costs
		"""
		print(self._column_costs, 'mean: ', sum(self._column_costs)/float(len(self._column_costs)))
		return self._column_costs

	def get_normalized_column_costs(self):
		"""
		Return the column costs normalized by number of reads in each column
		"""
		return self._normalized_column_costs


	def _compute_costs(self): 
		n = len(self._current_column)
		# first, get the clustering from the matrixcolumn
		id_to_names = defaultdict(list)
		for element in self._current_column:
			read = element[1]
			for var in read:
				if var.position == self._position:
					cluster = var.allele
					id_to_names[cluster].append(read.name)
		# printing
		print('computed column clustering:', self._positions)
		for k,v in id_to_names.items():
			print(k,v)

		# store for each true haplotype the computed cluster assignments
		individual_to_id = {'NA19240_HAP1':0, 'NA19240_HAP2':1, 'HG00514_HAP1':2, 'HG00514_HAP2':3}
		trueHT_to_computedHT = defaultdict(list)

		cluster_ids = list(id_to_names.keys())
		for i,cluster_id in enumerate(cluster_ids):
			for read_name in id_to_names[cluster_id]:
				readHT = read_name[-12:]
				trueHT_to_computedHT[individual_to_id[readHT]].append(i)

		# compute costs for flipping entries of same haplotype to the same allele
		HT_to_cost = defaultdict(lambda:[0]*self._ploidy)
		for trueHT,entries in trueHT_to_computedHT.items():
			for cluster in range(0,self._ploidy):
				# flip all entries to 'cluster'
				HT_to_cost[trueHT][cluster] += sum([0 if i==cluster else 1 for i in entries])

		# consider all possible assignments of 0,1,2,..,ploidy-1 to the partitions
		column_costs = []
		for a in itertools.permutations([i for i in range(0,self._ploidy)]):
			cost = 0
			for index,allele in enumerate(a):
				cost += HT_to_cost[index][allele]
			column_costs.append(cost)

		# printing
		print('true column clustering:', trueHT_to_computedHT)
		print('costs for current column:',min(column_costs))
		print('normalized costs for current column (cost/coverage):', min(column_costs)/float(n))
		self._column_costs.append(min(column_costs))
		self._normalized_column_costs.append(min(column_costs)/float(n))

