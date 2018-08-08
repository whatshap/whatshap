from whatshap.core import ReadSet, PedigreeDPTable, Pedigree, NumericSampleIds, PhredGenotypeLikelihoods
from whatshap.testhelpers import string_to_readset, brute_force_phase
from whatshap.phase import find_components
from whatshap.readsetpruning import ReadSetPruning
from whatshap.readsetpruning import ConflictSet

def generate_input(reads, weights=None):
	readset = string_to_readset(reads, weights)
	positions = readset.get_positions()
	components = find_components(positions, readset)
	return readset, positions, components

def check_conflict_set(column_clusters, expected_clustering):
	# construct ConflictSet
	unique_read_ids = set()
	for c in column_clusters:
		for e in c:
			for v in e:
				unique_read_ids.add(v)
	read_ids = list(unique_read_ids)
	conflict_set = ConflictSet(read_ids)
	for clustering in column_clusters:
		# compute mapping of elements to cluster_ids
		read_to_cluster = {}
		cluster_reads = []
		for i,cluster in enumerate(clustering):
			for read in cluster:
				read_to_cluster[read] = i
				cluster_reads.append(read)
		# add the conflicts/relationship to the ConflictSet
		n = len(cluster_reads)
		for id1 in range(n):
			for id2 in range(id1+1, n):
				read1 = cluster_reads[id1]
				read2 = cluster_reads[id2]
				if read_to_cluster[read1] == read_to_cluster[read2]:
					conflict_set.add_relationship(read1, read2)
				else:
					conflict_set.add_conflict(read1, read2)
	# check if clusters are computed correctly
	computed_clustering = sorted(conflict_set.get_clusters())
	print("expected clustering: ", expected_clustering)
	print("computed clustering: ", computed_clustering)
	assert computed_clustering == expected_clustering
				

def check_pruned_consensus_set(readset, expected_reads, expected_qualities):
	assert len(expected_reads) == len(expected_qualities)
	true_set = set([ (expected_reads[i],expected_qualities[i]) for i in range(len(expected_reads)) ])
	computed_set = set()
	
	for read in readset:
		# get allele sequence
		alleles = ''.join( str(var.allele) for var in read )
		qualities = ''.join( str(var.quality) for var in read )
		computed_set.add( (alleles,qualities) )

	print('computed reads: ', computed_set)
	print('true reads: ', true_set)
	assert(true_set == computed_set)


def test_cluster_merging1():
	column_clustering = [ [[1,2],[3]] , [[3,4], [5]]]
	expected_combined_clustering = [ [1,2], [3,4], [5] ]
	check_conflict_set(column_clustering, expected_combined_clustering)


def test_cluster_merging2():
	column_clustering = [ [[1,2]], [[2],[3]] ]
	expected_combined_clustering = [[1,2],[3]]
	check_conflict_set(column_clustering, expected_combined_clustering)


def test_cluster_merging3():
	column_clustering = [ [[1,2,5],[3,4]], [[1,3], [2,4]] ]
	expected_combined_clustering = [ [1,5], [2], [3], [4] ]
	check_conflict_set(column_clustering, expected_combined_clustering)


def test_cluster_merging4():
	column_clustering = [ [[1,2],[3]], [[3,4],[5]], [[3],[4,5]] ]
	expected_combined_clustering = [ [1,2], [3], [4], [5] ]
	check_conflict_set(column_clustering, expected_combined_clustering)


def test_cluster_merging5():
	column_clustering = [ [[1,2],[3]], [[3,4],[5]], [[3,5]] ]
	expected_combined_clustering = [[1,2],[3,4],[5]]
	check_conflict_set(column_clustering, expected_combined_clustering)


def test_cluster_merging6():
	column_clustering = [ [[1,2],[3]], [[3],[4,5]], [[1,4]] ]
	expected_combined_clustering = [ [1,2,4,5], [3] ]
	check_conflict_set(column_clustering, expected_combined_clustering)


def test_clustering1():
	reads = """
                11011011
                11111111
                00100100
                00101000
                00101000
                """
	weights = """
                11111211
                11211111
                11112211
                11111111
                11121111
                """
	readset, positions, components = generate_input(reads, weights)
	expected_clusters = [ ['Read 1', 'Read 2'], ['Read 3', 'Read 4', 'Read 5'] ]
	pruner = ReadSetPruning(readset, components, 2, 3, 6)
	computed_clusters = pruner.get_clusters()
	assert computed_clusters == expected_clusters

# TODO: disadvantage of the approach: Since only reads in a window are considered, and not all covering a region at the same time
#	the clustering misses some of the relationships between the reads. Here Read3 and Read8 are not combined, because they
#	are never within the same window (and no other reads in between relate to them)
def test_clustering2():
	reads = """
                101101111101101101110
                111101111101101111110111111
                001100001101010100000111
                111110011110111100110000011
                000000101001111000001110
                  1000000001111000001101111
                   110011111101100110010011
                        1100010100010111100
	        """
	readset, positions, components = generate_input(reads)
#	expected_clusters = [ ['Read 1', 'Read 2'], ['Read 3', 'Read 8'], ['Read 4', 'Read 7'], ['Read 5', 'Read 6'] ]
	expected_clusters = [ ['Read 1', 'Read 2'], ['Read 3'], ['Read 4', 'Read 7'], ['Read 5', 'Read 6'], ['Read 8'] ]
	pruner = ReadSetPruning(readset, components, 4, 5, 4)
	computed_clusters = pruner.get_clusters()
	assert computed_clusters == expected_clusters

def test_clustering3():
	reads = """
                111011100100111010100001000110101110
                110110100010111010100101000110001010
                111001111101010010010001110101100101
                111001101101010000000101110101000111
                010110001000111000110001010111110100
                000110000010111000110001011111010000
		"""

	readset, positions, components = generate_input(reads)
	expected_clusters = [ ['Read 1', 'Read 2'], ['Read 3', 'Read 4'], ['Read 5', 'Read 6']]
	pruner = ReadSetPruning(readset, components, 3, 5, 4)
	computed_clusters = pruner.get_clusters()
	assert computed_clusters == expected_clusters

def test_clustering4():
	reads = """
                11
                  00
                    10
		"""
	readset, positions, components = generate_input(reads)
	expected_clusters = [['Read 1'], ['Read 2'], ['Read 3']]
	pruner = ReadSetPruning(readset, components, 2, 2, 2)
	computed_clusters = pruner.get_clusters()
	assert sorted(computed_clusters) == sorted(expected_clusters)

def test_clustering5():
	reads = """
                11111
                 11110
                  10011
                       0100
                       01011
                        0000
		"""
	readset, positions, components = generate_input(reads)
	expected_clusters = [['Read 1', 'Read 2'], ['Read 3'], ['Read 4', 'Read 5'], ['Read 6']]
	pruner = ReadSetPruning(readset, components, 2, 3, 3)
	computed_clusters = pruner.get_clusters()
	print("computed clusters: ", computed_clusters)
	print("expected clusters: ", expected_clusters)
	assert sorted(computed_clusters) == sorted(expected_clusters)

# TODO: same issue as before: especially as ploidy increases, the number of reads per window needs to be large in order to not miss to many
#	relationships between reads.
def test_clustering6():
	reads = """
                1010001001101010001010100010001111111010100100
                1110100010001000101000000101000101100100011101000
                000010010101101100010001111000011000011000010000011
                000011111111100000100000011100001111111111011111001111
                101000100000000010001111100000010000100011100011111010
                101000000000000010001111101000000010000011100011111010
                010101010000000011100111111100011100001111110001011110
                010101010000000000000111111111111100001111110001111010
                000011000111000111001001111100001010101010001001001001
                001010100101010000000100101101101110111100011001011000
                001000100101010001100000101111101110111100011111011000
                   011110111100010000000000000001111111101111111001111
                    00010101111100010001111000000010111000000010011111
                    000010001010101000010101010100010101010100100111
                      000101000111001001111000101010101010001001001
                        1010101000100010100100111111101010010010001011
		"""
	readset, positions, components = generate_input(reads)
	expected_clusters = [['Read 1', 'Read 16'],['Read 14', 'Read 2'],['Read 13', 'Read 3'], ['Read 12', 'Read 4'],['Read 5', 'Read 6'],['Read 7', 'Read 8'],
				['Read 15', 'Read 9'], ['Read 10', 'Read 11'] ]
	pruner = ReadSetPruning(readset, components, 8, 20, 5)
	computed_clusters = pruner.get_clusters()
	print("computed: ", computed_clusters)
	print("expected: ", expected_clusters)
	assert computed_clusters == expected_clusters


# TODO: fails for number_of_clusters=2, since in the last window, it clusters the reads like this: {4,5}/{3}.
# 	Since similarity to read 3 is to different, the clusters are not combined.
#	To solve: #TODO better way to combine clusters that are similar
def test_clustering7():
	reads = """
                110110
                111111
                001001
                001010
                001110
		"""
	weights = """
                111112
                112111
                111122
                111111
                111211
		"""
	readset, positions, components = generate_input(reads, weights)
	expected_clusters = {	4:[ ['Read 1'], ['Read 2'], ['Read 3'], ['Read 4', 'Read 5'] ],
				3:[ ['Read 1', 'Read 2'], ['Read 3'], ['Read 4', 'Read 5'] ],
				2:[ ['Read 1', 'Read 2'], ['Read 3', 'Read 4', 'Read 5'] ],
				1:[ ['Read 1', 'Read 2', 'Read 3', 'Read 4', 'Read 5'] ]
				}

	for number_of_clusters in range(1,5):
		pruner = ReadSetPruning(readset, components, number_of_clusters, 4, 4)
		computed_clusters = pruner.get_clusters()
		print('computed clusters:', computed_clusters)
		print('expected clusters:', expected_clusters[number_of_clusters])
		assert computed_clusters == expected_clusters[number_of_clusters]

def test_clustering8():
	reads = """
                11
                1111
                  1111
                    11
		"""
	readset, positions, components = generate_input(reads)
	expected_clusters = [ ['Read 1', 'Read 2', 'Read 3', 'Read 4'] ]
	pruner = ReadSetPruning(readset, components, 2, 4, 2)
	computed_clusters = pruner.get_clusters()
	print("computed: ", computed_clusters)
	print("expected: ", expected_clusters)
	assert computed_clusters == expected_clusters


#def test_pruning_consensus1():
#	reads = """
#		110110
#		111111
#		001001
#		001010
#		001110
#		"""
#	weights = """
#		111112
#		112111
#		111122
#		111111
#		111211
#		"""
#	readset, positions, components = generate_input(reads, weights)
#	expected_reads = {
#			4:['110110','111111', '001001', '001110'],
#			3:['111110','001001','001110'],
#			2:['111110','001'],
#			1:['001110']
#			}
#
#	expected_qualities = {
#			4:['111112','112111','111122','222122'],
#			3:['221221','111122','222122'],
#			2:['221221','333'],
#			1:['114221']
#			}
#
#	for number_of_clusters in range(1,5):
#		pruner = ReadSetPruning(readset, components, number_of_clusters, 5, 3)
#		pruned_readset = pruner.get_pruned_readset()
#		print(pruned_readset)
#		check_pruned_consensus_set(pruned_readset, expected_reads[number_of_clusters], expected_qualities[number_of_clusters])
#
#
#def test_pruning_consensus2():
#	reads = """
#		11111
#		 11110
#		  10011
#		       0100
#		       01011
#		        0000
#		"""
#	readset, positions, components = generate_input(reads)
#	expected_reads = {
#			2:['111110', '10011', '0101', '0000'],
#			4:['111110', '10011', '0100', '01011', '0000']
#			}
#	expected_qualities = {
#			2:['122221','11111', '2221', '1111'],
#			4:['122221', '11111', '1111', '11111', '1111']
#			}
#
#	for number_of_clusters in [2,4]:
#		pruner = ReadSetPruning(readset, components, number_of_clusters, 3, 3)
#		pruned_readset = pruner.get_pruned_readset()
#		print(pruned_readset)
#		check_pruned_consensus_set(pruned_readset, expected_reads[number_of_clusters], expected_qualities[number_of_clusters])
#
#def test_pruning_consensus3():
#	reads = """
#		1100
#		00111111
#		    1111
#		"""
#	readset, positions, components = generate_input(reads)
#	pruner = ReadSetPruning(readset, components, 1, 2, 3)
#	pruned_readset = pruner.get_pruned_readset()
#	print(pruned_readset)
#
#	# since similarities between reads 2-1 and 3-1 are 0, these reads should not appear
#	# in the same clusters
#
#	expected_reads = ['1100','00111111']
#	expected_qualities = ['1111', '11112222']
#	check_pruned_consensus_set(pruned_readset, expected_reads, expected_qualities)
#
#def test_pruning_consensus4():
#	reads = """
#		1010
#		01011111
#		    0000
#		"""
#
#	readset, positions, components = generate_input(reads)
#	pruner = ReadSetPruning(readset, components, 1, 2, 3)
#	pruned_readset = pruner.get_pruned_readset()
#	print(pruned_readset)
#
#	expected_reads = ['1010', '01011111', '0000']
#	expected_qualities = ['1111', '11111111', '1111']
#	check_pruned_consensus_set(pruned_readset, expected_reads, expected_qualities)
#
# TODO in this case the consensus read is completely empty, since at each position
# both alleles have the same weight
# TODO similarly, what should happen if consensus read covers only a single variant
# (for same reason as before)
#def test_pruning_consensus5():
#	reads = """
#		0101
#		0111
#		1010
#		"""
#	weights = """
#		1121
#		1111
#		2212
#		"""
#	readset, positions, components = generate_input(reads, weights)
#	pruner = ReadSetPruning(readset, components, 1, 3, 3)
#	pruned_readset = pruner.get_pruned_readset()
#	print(pruned_readset)
#	# currently, consensus reads with less than 2 variants are not added
#	# to final set of reads
#	assert(len(pruned_readset) == 0)
#
#def test_pruning_consensus6():
#	reads = """
#		0101
#		0111
#		1010
#		"""
#	weights = """
#		1121
#		1111
#		2211
#		"""
#	readset, positions, components = generate_input(reads, weights)
#	pruner =  ReadSetPruning(readset, components, 1, 3, 3)
#	pruned_readset = pruner.get_pruned_readset()
#	print(pruned_readset)
#	# there should be no consensus since read would only contain a single variant
#	assert(len(pruned_readset) == 0)

