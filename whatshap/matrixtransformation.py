from collections import defaultdict
import sys
from .graph import ComponentFinder
from .core import Read, ReadSet, ClusterEditingSolver, DynamicSparseGraph
from .readscoring import score_global
import logging
import pysam
import itertools

logger = logging.getLogger(__name__)

class MatrixTransformation:
	"""
	Construct transformed matrix by computing clusterings for each variant positon.
	"""
	def __init__(self, reads, position_to_component, ploidy, min_overlap):
		"""
		reads -- ReadSet with reads to be considered
		position_to_component -- dict mapping variant positions to connected components
		ploidy -- the ploidy of the data
		min_overlap -- minimum required overlap of two reads
		"""
		# readset containing all reads of current connected component
		self._component_reads = None
		# set of reads covering current position
		self._current_column = None
		# cluster matrix aka transformed matrix
		self._cluster_matrix = ReadSet()
		# current mapping: read_name -> Read() object containing assignment to partitions
		self._readname_to_partitions = {}
		# number of clusters in each column
		self._number_of_clusters = []

		# map component_id -> read_ids
		connected_components = defaultdict(list)
		for i,read in enumerate(reads):
			position = read[0].position
			component_id = position_to_component[position]
			connected_components[component_id].append(i)
		# consider connected components separately
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
			num_clusters = []
			for j, position in enumerate(component_positions):
				self._current_column = []
				# get all reads that cover current position
				for i,read in enumerate(self._component_reads):
					if position in read:
						self._current_column.append( (i, read) )		
				if len(self._current_column) == 0:
					continue
				column = self._component_reads.subset([r[0] for r in self._current_column])

				# compute similarities for reads in column
				similarities = score_global(column, ploidy, min_overlap)
				
				# create read graph object
				graph = DynamicSparseGraph(len(column))
	
				# insert edges into read graph
				for (read1, read2) in similarities:
					graph.addEdge(read1, read2, similarities.get(read1, read2))
				
				# run cluster editing
				clusterediting = ClusterEditingSolver(graph, False)
				readpartitioning = clusterediting.run()
				num_clusters.append(len(readpartitioning))

				# store the result in final MEC matrix
				for c,cluster in enumerate(readpartitioning):
					for read_id in cluster:
						# get readname
						read_name = column[read_id].name
						quality = [10] * max(ploidy,len(readpartitioning))
						quality[c] = 0
						self._readname_to_partitions[read_name].add_variant(position, c, quality)
				self._number_of_clusters.append(max(ploidy,len(readpartitioning)))
				
			avg_clusters = sum(num_clusters)/len(component_positions)
			lt_k = sum([1 for cluster in num_clusters if cluster > ploidy])
			ht_k = sum([1 for cluster in num_clusters if cluster < ploidy])
			logger.debug("Transformed matrix of {} clustering with {} clusters on average ({} < k, {} > k)".format(len(num_clusters), avg_clusters, lt_k, ht_k))
			
			# add the computed reads to final ReadSet
			for read in self._readname_to_partitions.values():
				if len(read) > 0:
					self._cluster_matrix.add(read)
		# sort readsets
		self._cluster_matrix.sort()

	def get_transformed_matrix(self):
		"""
		Return the readset which represents the transformed matrix.
		"""
		for read in self._cluster_matrix:
			prev = -1
			for var in read:
				if var.position == prev:
					print('DOUBLE: ',read)
				prev = var.position
		return self._cluster_matrix

	def get_cluster_counts(self):
		"""
		Return the number of clusters per column in the cluster matrix.
		"""
		return self._number_of_clusters
