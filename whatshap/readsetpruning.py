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

class DotWriter:
	"""
	Writes a graph representation in .dot format.
	"""

	def __init__(self, filename):
		"""
		filename -- how to name the output file
		"""
		# filename
		self._filename = filename
		# edges of the graph represented as (node1, node2, label)
		self._edges = []

	def add_edge(self, node1, node2, label):
		"""
		Add an edge between node1 and node2 with given label.
		"""
		self._edges.append( (str(node1), str(node2), str(label)) )

	def write(self):
		start_nodes = set([i[0] for i in self._edges])
		end_nodes =  set([i[1] for i in self._edges]) - start_nodes
		all_nodes = start_nodes.union(end_nodes)
		level_nodes = [node for node in all_nodes if node.startswith('level')]
		list_end_nodes = ",".join(end_nodes)

		outfile = open(self._filename, 'w')
		outfile.write('graph graphname {\n')

		outfile.write('{ node [shape=plaintext]; \n ')
		for node in level_nodes:
			splitted = node.split(':')
			outfile.write(splitted[0] + ' [label=' + splitted[1] + '];')
		outfile.write( '\n' + "--".join(reversed([str(i+1) for i in range(len(level_nodes))])) + '\n')
		outfile.write('}')
		for edge in self._edges:
			start = edge[0].split(':')[0] if edge[0].startswith('level') else edge[0]
			end = edge[1].split(':')[0] if edge[1].startswith('level') else edge[1]
			outfile.write(start + '--' + end + ';\n')
		outfile.write('{rank = same;' + list_end_nodes + '}')
		for i in range(len(level_nodes)):
			outfile.write('{rank = same;' + str(i+1) + ', level' + str(i) + ' }\n')
		outfile.write('}')
		outfile.close()

class ReadSetPruning:
	"""
	Prune the ReadSet.
	
	Given all reads and positions to be used for phasing/genotyping,
	compute pairwise similarities and combine reads with high similarity
	into a single, consensus read.
	This is useful in order to reduce the state space.
	"""
	def __init__(self, reads, position_to_component, number_of_clusters, reads_per_window, variants_per_window, compute_consensus=False, write_dot = False):
		"""
		reads -- ReadSet with reads to be considered
		position_to_component -- dict mapping variant positions to connected components
		number_of_clusters -- how many steps should be performed for clustering
		reads_per_window -- max number of reads to consider per window
		variants_per_window -- min number of variants to consider per window
		compute_consensus -- if True, the clustered reads are combined to a consensus read
		"""
		# given parameters
		self._number_of_clusters = number_of_clusters
		self._reads_per_window = reads_per_window
		self._variants_per_window = variants_per_window
		self._compute_consensus = compute_consensus
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
		# clusters of current connected component
		self._clusters = []
		# conflict set of current connedted component
		self._conflict_set = None
		# pruned readset
		self._pruned_reads = ReadSet()
		# computed clusters as lists of read names (used for testing)
		self._readname_clusters = []
		# write column trees in dot file
		self._write_dot = write_dot
	
		# store overall mapping (readname->cluster)
		self._readname_to_cluster = {}
	
		# map component_id -> read_ids
		connected_components = defaultdict(list)
		for i,read in enumerate(reads):
			position = read[0].position
			component_id = position_to_component[position]
			connected_components[component_id].append(i)
		for component_id, component in connected_components.items():
			# get reads belonging to this component
			self._component_reads = reads.subset(component)
			# get all positions in this component
			component_positions = self._component_reads.get_positions()
			# create a ConflictSet for the reads in this component
			self._conflict_set = ConflictSet([i for i in range(len(self._component_reads))])

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
					k += 1

				# compute similarities
				self._compute_similarities()
#				print('similarities for window: ', self._window, self._similarities)
				# cluster based on similarities
				self._compute_clusters()
				self._window += 1
				print('current overall clustering:')
				print([[self._component_reads[j].name for j in i ] for i in self._conflict_set.get_clusters()])

#				if (i + self._reads_per_window) == len(self._component_reads):
#					break

			# get the computed read clusters
			clusters = self._conflict_set.get_clusters()
			for c in clusters:
				# construct readset
				self._clusters.append(self._component_reads.subset(c))

			# store for each read, the cluster it belongs to
			for i, c in enumerate(self._clusters):
				cluster_id = 'connected_component_' + str(component_id) + '_cluster_' + str(i)
				for read in c:
					self._readname_to_cluster[read.name] = cluster_id

			print('readname_to_cluster:', self._readname_to_cluster)

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

		dotwriter = DotWriter('window_' + str(self._window) + '.dot')
		cluster_to_name = defaultdict(lambda : None)

		# perform k steps
#		k = n - self._number_of_clusters
		k = n
		logger.debug('number of iterations: %d', k)
		
		while (k > 0):
#		for i in range(k):
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
			c1 = clusters.find(max_column)
			c2 = clusters.find(max_row)
			clusters.merge(max_column, max_row)

			## write dot file
			if self._write_dot:
				node_label = 'level' + str(i)+':'+str(max_value)
				dotwriter.add_edge(node_label,cluster_to_name[c1] if cluster_to_name[c1] is not None else self._current_column[c1][1].name, max_value)
				dotwriter.add_edge(node_label,cluster_to_name[c2] if cluster_to_name[c2] is not None else self._current_column[c2][1].name, max_value)
				cluster_to_name[clusters.find(max_column)] = node_label

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

#		id_to_index = defaultdict(list)
#		# collect all read ids beloning to same cluster
#		for i in range(n):
#			id_to_index[clusters.find(i)].append(self._current_column[i][0])
#
#		# compute consensus strings
#		consensus_strings = {}
#		for cluster_id, cluster in id_to_index.items():
#			# get reads of this cluster
#			cluster_reads = self._component_reads.subset(cluster)
#			# compute consensus based on current positions # TODO
#			consensus_strings[cluster_id] = self._consensus_string(cluster_reads)
#		print('consensus strings: ', consensus_strings)
#		# combine clusters that have identical consensus strings
#		cluster_ids = list(id_to_index.keys())
#		for i in range(len(cluster_ids)):
#			for j in range(i+1, len(cluster_ids)):
#				cluster_id_i = cluster_ids[i]
#				cluster_id_j = cluster_ids[j]
#				if consensus_strings[cluster_id_i] == consensus_strings[cluster_id_j]:
#					# merge clusters
#					clusters.merge(id_to_index[cluster_id_i][0], id_to_index[cluster_id_j][0])

		# printing
		id_to_names = defaultdict(list)
		for i in range(n):
			id_to_names[clusters.find(i)].append(self._current_column[i][1].name)

		print('current column clustering:', self._positions)
		for k,v in id_to_names.items():
			print(k,v)

		# write dot file
		if self._write_dot:
			dotwriter.write()

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
					if variant.allele == 0:
						pos_to_allele[variant.position] += variant.quality
					elif variant.allele == 1:
						pos_to_allele[variant.position] -= variant.quality
			# add variants to consensus read
			for pos in pos_to_allele.keys():
				quality = pos_to_allele[pos]
				if quality != 0:
					consensus_allele = 0 if quality > 0 else 1
					consensus_read.add_variant(pos, consensus_allele, abs(quality))
			consensus_read.sort()
			if len(consensus_read) > 1:
				self._pruned_reads.add(consensus_read)

#	def _consensus_string(self, readset):
#		"""
#		Compute consensus for currently considered positions for the reads contained
#		in the given readset.
#		"""
#		pos_to_allele = defaultdict(int)
#		consensus_string = ''
#		for read in readset:
#			for variant in read:
#				if variant.position not in self._positions:
#					continue
#				for i in range(len(variant.allele)):
#					if variant.allele[i] == 0:
#						pos_to_allele[variant.position] += variant.quality[i]
#					elif variant.allele[i] == 1:
#						pos_to_allele[variant.position] -= variant.quality[i]
#		for pos in pos_to_allele.keys():
#			quality = pos_to_allele[pos]
#			if quality != 0:
#				consensus_allele = 0 if quality > 0 else 1
#				consensus_string += str(consensus_allele)
#			else:
#				consensus_string += '3'
#		return consensus_string


	def tag_reads(self, bam_file, output_filename):
		"""
		In the given BAM file, tag the reads according to the clusters they have been assigned to.
		"""
		# TODO implement this
		pass

		bam_reader = pysam.AlignmentFile(bam_file)
		bam_writer = (output_filename, 'wb')
		for alignment in bam_reader:
			cluster_id = self._readname_to_cluster[alignment.query_name]
			alignment.set_tag('XC', value=cluster_id)
			bam_writer.write(alignment)
		bam_reader.close()
		bam_writer.close()
