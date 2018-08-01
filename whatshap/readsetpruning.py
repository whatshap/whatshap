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

class ConflictSet:
	"""
	For a set of reads, keep track of which reads are in conflict
	and must not end up in the same cluster.
	"""
	def __init__(self, read_ids):
		self._read_ids = sorted(read_ids)
		self._conflicts = defaultdict(int)

	def add_conflict(self, read1, read2):
		"""
		Add a conflict between the given reads.
		"""
		assert (read1 in self._read_ids) and (read2 in self._read_ids)
		if read1 < read2:
			self._conflicts[(read1,read2)] = 1
		else:
			self._conflicts[(read2,read1)] = 1

	def add_relationship(self, read1, read2):
		"""
		Mark the relationship between the reads, i.e. if they are in conflict,
		nothing changes, if not, this indicates that the reads are currently in the
		same cluster.
		"""
		assert (read1 in self._read_ids) and (read2 in self._read_ids)
		if read1 < read2:
			if not self.in_conflict(read1,read2):
				self._conflicts[(read1,read2)] = -1
		else:
			if not self.in_conflict(read1,read2):
				self._conflicts[(read2,read1)] = -1

	def in_conflict(self, read1, read2):
		"""
		Check if the given reads are in conflict.
		"""
		assert (read1 in self._read_ids) and (read2 in self._read_ids)
		if read1 < read2:
			return self._conflicts[(read1,read2)] == 1
		else:
			return self._conflicts[(read2,read1)] == 1

	def in_same_cluster(self, read1, read2):
		"""
		Check if the given reads are in the same cluster.
		"""
		assert (read1 in self._read_ids) and (read2 in self._read_ids)
		if read1 < read2:
			return self._conflicts[(read1,read2)] == -1
		else:
			return self._conflicts[(read2,read1)] == -1

#	def get_clusters(self):
#		"""
#		Return lists of objects that are not in conflict.
#		"""
#		result = []
#		used_reads = []
#		n = len(self._read_ids)
#		for i in range(n):
#			if i in used_reads:
#				continue
#			used_reads.append(i)
#			cluster = [self._read_ids[i]]
#			for j in range(i+1, n):
#				if j in used_reads:
#					continue
#				if self.in_same_cluster(self._read_ids[i], self._read_ids[j]):
#					cluster.append(self._read_ids[j])
#					used_reads.append(j)
#			result.append(cluster)
#		print(self._conflicts)
#		return result

	def get_clusters(self):
		result = []
		used_reads = []
		n = len(self._read_ids)
		for i in range(n):
			if i in used_reads:
				continue
			used_reads.append(i)
			next_used_reads, cluster = self.get_cluster(i, used_reads, [self._read_ids[i]], n)
			used_reads = next_used_reads
			result.append(cluster)
		return sorted(result)

	def get_cluster(self, i, used_reads, cluster, n):
		for j in range(i+1,n):
			if j in used_reads:
				continue
			if self.in_same_cluster(self._read_ids[i], self._read_ids[j]):
				conflict = False
				for c in cluster:
					if c == self._read_ids[j]:
						continue
					if self.in_conflict(self._read_ids[j], c):
						conflict = True
						break
				if not conflict:
					cluster.append(self._read_ids[j])
					used_reads.append(j)
					next_used_reads, next_cluster = self.get_cluster(j, used_reads, cluster, n)
					cluster = next_cluster
					used_reads = next_used_reads
		return used_reads, cluster

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
		number_of_clusters -- how many steps should be performed for clustering
		compute_consensus -- if True, the clustered reads are combined to a consensus read
		"""

		self._number_of_clusters = number_of_clusters
		self._compute_consensus = compute_consensus
		# currently considered position
		self._position = None
		# set of reads covering current position
		self._current_column = None
		# similarities of current position
		self._similarities = None
		# clusters of current connected component
		self._clusters = []
		# conflict set of current connedted component
		self._conflict_set = None
		# pruned readset
		self._pruned_reads = ReadSet()
		# computed clusters as lists of read names (used for testing)
		self._readname_clusters = []
		
		# map component_id -> read_ids
		connected_components = defaultdict(list)
		for i,read in enumerate(reads):
			position = read[0].position
			component_id = position_to_component[position]
			connected_components[component_id].append(i)
		for component_id, component in connected_components.items():
			# get reads belonging to this component
			component_reads = reads.subset(component)
			# get all positions in this component
			component_positions = component_reads.get_positions()
			# create a ConflictSet for the reads in this component
			self._conflict_set = ConflictSet([i for i in range(len(component_reads))])

			previous_reads = None
			for pos in component_positions:
				# find all reads covering this position
				self._current_column = [(i,read) for i,read in enumerate(component_reads) if pos in read]
				read_ids = [e[0] for e in self._current_column]
				# if column is covered by same reads as previous, skip it
				if previous_reads is not None:
					if previous_reads == read_ids:
						continue
				previous_reads = read_ids
				# compute similarities
				self._compute_similarities()
				# cluster based on similarities
				self._compute_clusters()
#				print('current_clustering', self._conflict_set.get_clusters())

			# get the computed read clusters
			clusters = self._conflict_set.get_clusters()
			for c in clusters:
				# construct readset
				self._clusters.append(component_reads.subset(c))
			if self._compute_consensus:
				self._consensus_reads()
			else:
				self._combine_reads()
			# store string representation of clusters
			for readset in self._clusters:
				self._readname_clusters.append(sorted([read.name for read in readset]))
			self._clusters = []

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
		n = len(self._current_column)

		# precompute lists of variants supported by each read
		variants_per_read = []
		for i in range(n):
			variants = { variant.position:variant.allele for variant in self._current_column[i][1] }
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
		id_to_names = defaultdict(list)
		for i in range(n):
			id_to_names[clusters.find(i)].append(self._current_column[i][0])


#		print('current column clustering:')
#		for k,v in id_to_names.items():
#			print(k,v)

		# based on the clustering, update the ConflictSet
		for i in range(n):
			for j in range(i+1,n):
				if clusters.find(i) is not clusters.find(j):
					self._conflict_set.add_conflict(self._current_column[i][0],self._current_column[j][0])
				else:
					self._conflict_set.add_relationship(self._current_column[i][0],self._current_column[j][0])

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
