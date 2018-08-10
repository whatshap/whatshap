from collections import defaultdict
import sys
from .graph import ComponentFinder
from .core import Read, ReadSet
import logging
import pysam

logger = logging.getLogger(__name__)

def compute_similarity(variants1, variants2):
	"""
	Compute a similarity score for a given pair of reads.
	"""
	# determine which positions overlap
	positions1 = set(variants1.keys())
	positions2 = set(variants2.keys())
	overlap = positions1.intersection(positions2)
	if len(overlap) < 2:
		return 0.0
	shared_alleles = 0
	for pos in overlap:
		if variants1[pos] == variants2[pos]:
			shared_alleles += 1
	return shared_alleles / len(overlap)


class ReadSetPruning:
	"""
	Prune the ReadSet.
	
	Given all reads and positions to be used for phasing/genotyping,
	compute pairwise similarities and combine reads with high similarity
	into a single, consensus read.
	This is useful in order to reduce the state space.
	"""
	def __init__(self, reads, position_to_component, number_of_clusters, reads_per_window, variants_per_window):
		"""
		reads -- ReadSet with reads to be considered
		position_to_component -- dict mapping variant positions to connected components
		number_of_clusters -- how many steps should be performed for clustering
		reads_per_window -- max number of reads to consider per window
		variants_per_window -- min number of variants to consider per window
		"""
		# given parameters
		self._number_of_clusters = number_of_clusters
		self._reads_per_window = reads_per_window
		self._variants_per_window = variants_per_window
		# currently considered positions
		self._positions = None
		# currently considered window
		self._window = 0
		# readset containing all reads of current connected component
		self._component_reads = None
		# set of reads covering current position
		self._current_column = None
		# similarities of current position
		self._similarities = None
		# cluster matrix
		self._cluster_matrix = ReadSet()
		# readset with used positions
		self._allele_matrix = ReadSet()
		# current mapping: read_name -> Read() object containing assignment to partitions
		self._readname_to_partitions = {}
		# current mapping: read name -> alleles used for comparison
		self._readname_to_positions = defaultdict(set)

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
			# create map containing positions used for each read
			self._readname_to_positions = defaultdict(set)
			for read in self._component_reads:
				self._readname_to_partitions[read.name] = Read(read.name, read.mapqs[0], read.source_id, read.sample_id)

			# get all positions in this component
			component_positions = self._component_reads.get_positions()
			# was the last read in any of the windows
			last_used = False
			# consider windows of reads
			for i, read in enumerate(self._component_reads):
				# add read to current column
				self._current_column = [(i,read)]
				# get variant positions covered by this read
				self._positions = set([var.position for var in read])
				# iterate over next k reads
				k = 1
				while ( k < self._reads_per_window) and ( i+k < len(self._component_reads)):
					# get intersection of variant positions
					current_read = self._component_reads[i+k]
					var_positions = set([var.position for var in current_read])
					intersection = self._positions.intersection(var_positions)
					# if intersection sufficiently large, consider read
					if len(intersection) >= self._variants_per_window:
						self._current_column.append( (i+k, current_read) )
						self._positions = intersection
						# if last read was used, don't look at next windows
						if (i+k) == len(self._component_reads) - 1:
							last_used = True
					else:
						break
 
					k += 1

				# if less than 2 reads are in the same window, skip it
				if len(self._current_column) < 2:
					continue

				# for each read in the window, remember used alleles
				for (i,read) in self._current_column:
					self._readname_to_positions[read.name] = self._readname_to_positions[read.name].union(self._positions)

				# compute similarities
				self._compute_similarities()
				# cluster based on similarities
				self._compute_clusters()
				self._window += 1
				if last_used:
					break

			# add the computed reads to final ReadSet
			for read in self._readname_to_partitions.values():
				if len(read) > 0:
					self._cluster_matrix.add(read)

			# construct allele matrix containing only positions used for comparing reads during clustering
			for read in self._component_reads:
				new_read = Read(read.name, read.mapqs[0], read.source_id, read.sample_id)
				for var in read:
					if var.position in self._readname_to_positions[read.name]:
						new_read.add_variant(var.position, var.allele, var.quality)
				if len(read) > 2:
					self._allele_matrix.add(new_read)
			
		# sort readsets
		self._cluster_matrix.sort()
		self._allele_matrix.sort()


	def get_cluster_matrix(self):
		"""
		Return the readset which represents the cluster matrix.
		"""
		return self._cluster_matrix

	def get_allele_matrix(self):
		"""
		Return the readset containing those positions for each read
		which have been used during clustering.
		"""
		return self._allele_matrix

	def _compute_similarities(self):
		"""
		Compute pairwise read similarities for a component.
		"""
		n = len(self._current_column)

		# precompute lists of variants supported by each read
		variants_per_read = []
		for i in range(n):
			variants = { variant.position:variant.allele for variant in self._current_column[i][1] if variant.position in self._positions}
			variants_per_read.append(variants)
		self._similarities = [ [ -1 for x in range(n) ] for y in range(n) ]
		for i in range(n):
			for j in range(i+1,n):
				self._similarities[i][j] = round(compute_similarity(variants_per_read[i], variants_per_read[j]),1)
				self._similarities[j][i] = self._similarities[i][j]

	def _compute_clusters(self):
		"""
		Perform hierarchical clustering to find groups of similar reads
		"""
		n = len(self._current_column)
		clusters = ComponentFinder([i for i in range(n)])

		# perform k steps
		k = n
		logger.debug('number of iterations: %d', k)
		
		while (k > 0):
			# determine which clusters to merge
			max_value = -float('inf')
			max_column = -1
			max_row = -1
			for a in range(n):
				for b in range(a+1, n):
					if self._similarities[a][b] > max_value:
						max_value = self._similarities[a][b]
						max_column = a
						max_row = b
			if max_value <= 0:
				# reads are not overlapping and should not be merged
				break
			if (k <= self._number_of_clusters):
				# only proceed if similarity is high enough
				# TODO: proper theshold
				if max_value < 0.8:
					break

			# combine clusters
			logger.debug('iteration %d:', n-k)
			logger.debug('merge read  %s and %s (similarity: %f)', self._current_column[max_column][1].name, self._current_column[max_row][1].name, max_value)
			clusters.merge(max_column, max_row)

			# recompute the similarities based on average linkage
			for j in range(n):
				if (j != max_row) and (j != max_column):
					average_dist = (self._similarities[max_row][j] + self._similarities[max_column][j]) / 2
					self._similarities[max_column][j] = average_dist
					self._similarities[j][max_column] = average_dist
			for j in range(n):
				self._similarities[j][max_row] = -float('inf')
				self._similarities[max_row][j] = -float('inf')
			k -= 1

		# printing
		id_to_names = defaultdict(list)
		for i in range(n):
			id_to_names[clusters.find(i)].append(self._current_column[i][1].name)

		print('current column clustering:', self._positions)
		for k,v in id_to_names.items():
			print(k,v)

		# store the clustering in MEC matrix
		cluster_ids = list(id_to_names.keys())
		for i,cluster_id in enumerate(cluster_ids):
			for read_name in id_to_names[cluster_id]:
				self._readname_to_partitions[read_name].add_variant(self._window, i, 10*len(self._positions))

