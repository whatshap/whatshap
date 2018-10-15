from collections import defaultdict
import sys
from .graph import ComponentFinder
from .core import Read, ReadSet, CoreAlgorithm, LightCompleteGraph
from .readscoring import score
import logging
import pysam
import itertools

logger = logging.getLogger(__name__)

def compute_similarity(variants1, variants2):
	"""
	Compute a similarity score for a given pair of reads.
	"""
	# determine which positions overlap
	positions1 = set(variants1.keys())
	positions2 = set(variants2.keys())
	overlap = positions1.intersection(positions2)
	if len(overlap) == 0:
		return 0.0
	shared_alleles = 0
	for pos in overlap:
		if variants1[pos] == variants2[pos]:
			shared_alleles += 1
	return shared_alleles / len(overlap)


class MatrixTransformation:
	"""
	Construct transformed matrix by computing clusterings for each variant positon.
	"""
	def __init__(self, reads, position_to_component, ploidy, errorrate, min_overlap):
		"""
		reads -- ReadSet with reads to be considered
		position_to_component -- dict mapping variant positions to connected components
		"""
		# readset containing all reads of current connected component
		self._component_reads = None
		# set of reads covering current position
		self._current_column = None
		# cluster matrix
		self._cluster_matrix = ReadSet()
		# current mapping: read_name -> Read() object containing assignment to partitions
		self._readname_to_partitions = {}
		# current mapping: read name -> alleles used for comparison
		self._readname_to_positions = defaultdict(set)
		# column costs
		self._column_costs = []
		self._number_of_clusters = []

		# map component_id -> read_ids
		connected_components = defaultdict(list)
		for i,read in enumerate(reads):
			position = read[0].position
			component_id = position_to_component[position]
			connected_components[component_id].append(i)
		for component_id, component in connected_components.items():
			# get reads belonging to this component
			self._component_reads = reads.subset(component)

			# create new Read objects that hold the cluster assignments
			self._readname_to_partitions = {}
			for read in self._component_reads:
				self._readname_to_partitions[read.name] = Read(read.name, read.mapqs[0], read.source_id, read.sample_id)
			# current column
			self._current_column = []
			# get all positions in this component
			component_positions = self._component_reads.get_positions()
			for j, position in enumerate(component_positions):
				self._current_column = []
				# get all reads that cover current position
				for i,read in enumerate(self._component_reads):
					if position in read:
						self._current_column.append( (i, read) )		
				if len(self._current_column) == 0:
					continue
				column = self._component_reads.subset([r[0] for r in self._current_column])
				print('reads covering position ', position, [r.name for r in column])

				# compute similarities for reads in column
				similarities = score(column, ploidy, errorrate, min_overlap)
#				print(similarities)

				# create read graph object
				graph = LightCompleteGraph(len(column),True)
	
				# insert edges into read graph
				n_reads = len(column)
				for id1 in range(n_reads):
					for id2 in range(id1+1, n_reads):
						graph.setWeight(id1, id2, similarities[id1][id2 - id1 - 1])

				# run cluster editing
				clusterediting = CoreAlgorithm(graph)	
				readpartitioning = clusterediting.run()
				print(readpartitioning)

				# store the result in final MEC matrix
				for c,cluster in enumerate(readpartitioning):
					for read_id in cluster:
						# get readname
						read_name = column[read_id].name
						quality = [10] * len(readpartitioning)
						quality[c] = 0
						print(read_name)
						self._readname_to_partitions[read_name].add_variant(position, c, quality)
				self._number_of_clusters.append(len(readpartitioning))

			# add the computed reads to final ReadSet
			for read in self._readname_to_partitions.values():
				if len(read) > 0:
					self._cluster_matrix.add(read)
		# sort readsets
		self._cluster_matrix.sort()
		print('')

	def get_transformed_matrix(self):
		"""
		Return the readset which represents the cluster matrix.
		"""
		return self._cluster_matrix

	def get_cluster_counts(self):
		"""
		Return the number of clusters per column in the cluster matrix.
		"""
		return self._number_of_clusters
