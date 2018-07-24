from collections import defaultdict
import sys
from .graph import ComponentFinder
from .core import Read, ReadSet
import logging

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
	def __init__(self, reads, position_to_component, number_of_clusters, compute_consensus=False):
		"""
		reads -- ReadSet with reads to be considered
		position_to_component -- dict mapping variant positions to connected components
		scaling_factor -- by which factor should a read component be reduced
		compute_consensus -- if True, the clustered reads are combined to a consensus read
		"""
		self._number_of_clusters = number_of_clusters
		self._compute_consensus = compute_consensus
		# connected components of read-overlap graph
		self._connected_components = defaultdict(list)
		# which component is currently considered
		self._component = None
		# similarities of current component
		self._similarities = None
		# clusters of current component
		self._clusters = None
		# pruned readset
		self._pruned_reads = ReadSet()
		self._readname_clusters = []

		# map index of component to reads it contains
		for read in reads:
			position = read[0].position
			component_id = position_to_component[position]
			self._connected_components[component_id].append(read)

		# prune the original set of reads
		for component_id in self._connected_components.keys():
			self._component = component_id
			self._compute_similarities()
			self._compute_clusters()
			print('computed clusters:')
			for readset in self._clusters:
				self._readname_clusters.append(sorted([read.name for read in readset]))
				print([read.name for read in readset])
			if self._compute_consensus:
				self._consensus_reads()
			else:
				self._combine_reads()

	# this just returns the computed read clusters as lists of read names
	# function mainly used for testing
	def get_clusters(self):
		return self._readname_clusters

	def get_pruned_readset(self):
		"""
		Return the new readset.
		"""
		self._pruned_reads.sort()
		return self._pruned_reads

	def _compute_similarities(self):
		"""
		Compute pairwise read similarities for a component.
		"""
		n = len(self._connected_components[self._component])

		# precompute lists of variants supported by each read
		variants_per_read = []
		for i in range(n):
			variants = { variant.position:variant.allele for variant in self._connected_components[self._component][i] }
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
		n = len(self._connected_components[self._component])
		clusters = ComponentFinder([i for i in range(n)])

		# perform k steps
#		k = int(round(n*float(self._scaling_factor),1))
		k = n - self._number_of_clusters
		logger.debug('number of iterations: %d', k)
		for i in range(k):
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
			
			# combine clusters
			logger.debug('iteration %d:', i)
			logger.debug('merge read  %d and %d (similarity: %d)', max_column, max_row, max_value)
#			for j in range(n):
#				print(self._similarities[j])
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

		id_to_cluster = defaultdict(list)
		for i in range(n):
			id_to_cluster[clusters.find(i)].append(i)

		# represent each cluster as a readset
		self._clusters = []
		for c in id_to_cluster:
			readset = ReadSet()
			for read in id_to_cluster[c]:
				readset.add(self._connected_components[self._component][read])
			readset.sort()
			self._clusters.append(readset)

	# TODO what happens if consensus read is empty or contains only one positon?
	# currently, such clusters are ignored (no consensus read is added to the final
	# readset, and the original reads are also not added)
	def _consensus_reads(self):
		"""
		Compute consensus reads from the cluters.
		"""
		# get reads belonging to a cluster
		for cluster in self._clusters:
			# count alleles at each position
			pos_to_allele = defaultdict(int)
			# the new consensus read
			consensus_read = Read(cluster[0].name, cluster[0].mapqs[0], cluster[0].source_id,
					cluster[0].sample_id, cluster[0].reference_start, cluster[0].BX_tag)
			# if read carries 0, add quality, o.w. substract
			for read in cluster:
				for variant in read:
					assert len(variant.allele) == len(variant.quality)
					for i in range(len(variant.allele)):
						if variant.allele[i] == 0:
							pos_to_allele[variant.position] += variant.quality[i]
						elif variant.allele[i] == 1:
							pos_to_allele[variant.position] -= variant.quality[i]
			# add variants to consensus read
			for pos in pos_to_allele.keys():
				quality = pos_to_allele[pos]
				if quality != 0:
					consensus_allele = 0 if quality > 0 else 1
					consensus_read.add_variant(pos, [consensus_allele], [abs(quality)])
			consensus_read.sort()
			if len(consensus_read) > 1:
				self._pruned_reads.add(consensus_read)

	# instead of computing consensus reads, group the clustered ones within the readset to keep
	# all alleles and qualities
	def _combine_reads(self):
		"""
		Combine the reads that were clustered together.
		"""
		for cluster in self._clusters:
			# Read object containing infomation of all reads of the cluster
			# TODO what to do with BX tag?
			combined_read = Read(cluster[0].name, cluster[0].mapqs[0], cluster[0].source_id, 
					cluster[0].sample_id, cluster[0].reference_start, cluster[0].BX_tag)
			
			# get alleles and qualities of reads
			pos_to_allele = defaultdict(list)
			pos_to_quality = defaultdict(list)

			for read in cluster:
				for variant in read:
					pos = variant.position
					pos_to_allele[pos].extend(variant.allele)
					pos_to_quality[pos].extend(variant.quality)
			# create a combined read
			for pos in pos_to_allele.keys():
				combined_read.add_variant(pos, pos_to_allele[pos], pos_to_quality[pos])
			combined_read.sort()
			if len(combined_read) > 1:
				self._pruned_reads.add(combined_read)
