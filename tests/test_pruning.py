from whatshap.core import ReadSet, PedigreeDPTable, Pedigree, NumericSampleIds, GenotypeLikelihoods, Genotype
from whatshap.testhelpers import string_to_readset, brute_force_phase
from whatshap.phase import find_components
from whatshap.readsetpruning import ReadSetPruning

from collections import defaultdict

def generate_input(reads, weights=None):
	readset = string_to_readset(reads, weights)
	positions = readset.get_positions()
	components = find_components(positions, readset)
	return readset, positions, components

def check_window_clustering(cluster_matrix, expected_clusters):
	"""
	Given cluster readset and list of expected
	column wise clusters. Compare the clusterings.
	"""

	position_to_clusters = defaultdict(lambda: defaultdict(list))
	for read in cluster_matrix:
		for var in read:
			position_to_clusters[var.position][var.allele].append(read.name)
	matrix_positions = list(position_to_clusters.keys())
	expected_positions = list(expected_clusters.keys())
	assert sorted(matrix_positions) == sorted(expected_positions)

	all_clusters = []
	all_expected_clusters = []
	for position in expected_positions:
		all_clusters.append(sorted([ sorted(cluster) for cluster in position_to_clusters[position].values()]))
		all_expected_clusters.append(sorted(expected_clusters[position]))
	assert sorted(all_clusters) == sorted(all_expected_clusters)

def solve_MEC(cluster_matrix, ploidy):
	"""
	finds the best consensus clustering solving MEC.
	"""
	numeric_sample_ids = NumericSampleIds()
	pedigree = Pedigree(numeric_sample_ids, ploidy)
	windows = cluster_matrix.get_positions()
	# allowed genotypes (require different allele for each partition, e.g. gt 0/1/2/3 for ploidy=4)
	genotypes = [Genotype([i for i in range(0,ploidy)])] * len(windows)
	pedigree.add_individual('0', genotypes, [GenotypeLikelihoods(ploidy, ploidy,[])]*len(windows))
	# n_alleles per window == ploidy
	allele_counts = [ploidy] * len(windows)
	dp_table = PedigreeDPTable(cluster_matrix, [1]*len(windows), pedigree, ploidy, distrust_genotypes=False, allele_counts=allele_counts)
	result = []
	for i in range(ploidy):
		result.append([])
	optimal_partitioning = dp_table.get_optimal_partitioning()
	for i,partition in enumerate(optimal_partitioning):
		result[partition].append(cluster_matrix[i].name)
	for i in range(len(result)):
		result[i] = sorted(result[i])
	return sorted(result), optimal_partitioning


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
	pruner = ReadSetPruning(readset, components, 2, 3, 6)
	cluster_matrix = pruner.get_cluster_matrix()
	allele_matrix = pruner.get_allele_matrix()
	print('cluster matrix: ', cluster_matrix)
	print('allele matrix: ', allele_matrix)

	expected_clusters = {
		0: [ ['Read 1', 'Read 2'], ['Read 3'] ],
		1: [ ['Read 2'], ['Read 3', 'Read 4'] ],
		2: [ ['Read 3', 'Read 4', 'Read 5'] ]
	}

	check_window_clustering(cluster_matrix, expected_clusters)
	# compute consensus clustering
	consensus_clustering, optimal_partitioning = solve_MEC(cluster_matrix, 2)
	expected_overall_clustering = sorted([ ['Read 1', 'Read 2'], ['Read 3', 'Read 4', 'Read 5'] ])
	print('expected clusters: ', expected_overall_clustering)
	print('computed clusters: ', sorted(consensus_clustering))
	assert sorted(consensus_clustering) == expected_overall_clustering

# TODO: disadvantage of the approach: Since only reads in a window are considered, and not all covering a region at the same time
#	the clustering misses some of the relationships between the reads. Here, if window size = 5, Read3 and Read8 are not combined, because they
#	are never within the same window (and no other reads in between relate to them). They are either put together with R1,2 or R3.
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
	pruner = ReadSetPruning(readset, components, 4, 6, 4)
	cluster_matrix = pruner.get_cluster_matrix()
	allele_matrix = pruner.get_allele_matrix()
	print('cluster matrix', cluster_matrix)
	print('allele_matrix', allele_matrix)

	expected_clusters = {
		0: [ ['Read 1', 'Read 2'], ['Read 3'], ['Read 4'], ['Read 5','Read 6'] ],
		1: [ ['Read 2'], ['Read 3'], ['Read 4','Read 7'], ['Read 5', 'Read 6'] ],
		2: [ ['Read 3','Read 8'], ['Read 4', 'Read 7'], ['Read 5', 'Read 6'] ]
	}

	check_window_clustering(cluster_matrix, expected_clusters)
	# compute consensus clustering
	consensus_clustering, optimal_partitioning = solve_MEC(cluster_matrix, 4)
	expected_overall_clustering = sorted([ ['Read 1', 'Read 2'], ['Read 3', 'Read 8'], ['Read 4', 'Read 7'], ['Read 5', 'Read 6'] ])
	print('expected clusters: ', expected_overall_clustering)
	print('computed clusters: ', sorted(consensus_clustering))
	assert sorted(consensus_clustering) == expected_overall_clustering

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
	pruner = ReadSetPruning(readset, components, 3, 5, 4)
	cluster_matrix = pruner.get_cluster_matrix()
	allele_matrix = pruner.get_allele_matrix()
	print('cluster matrix:', cluster_matrix)
	print('allele matrix:', allele_matrix)

	expected_clusters = {
		0: [ ['Read 1', 'Read 2'], ['Read 3', 'Read 4'], ['Read 5'] ],
		1: [ ['Read 2'], ['Read 3', 'Read 4'], ['Read 5', 'Read 6'] ]
	}

	check_window_clustering(cluster_matrix, expected_clusters)
	# get consensus clustering
	consensus_clustering, optimal_partitioning = solve_MEC(cluster_matrix, 3)
	expected_overall_clustering = sorted([ ['Read 1', 'Read 2'], ['Read 3', 'Read 4'], ['Read 5', 'Read 6'] ])
	print('expected clusters: ', expected_overall_clustering)
	print('computed clusters: ', sorted(consensus_clustering))
	assert sorted(consensus_clustering) == expected_overall_clustering


def test_clustering4():
	reads = """
                11
                  00
                    10
		"""
	# each window contains only one read, so they will be skipped
	readset, positions, components = generate_input(reads)
	pruner = ReadSetPruning(readset, components, 2, 2, 2)
	cluster_matrix = pruner.get_cluster_matrix()
	allele_matrix = pruner.get_allele_matrix()
	print('cluster_matrix: ', cluster_matrix)
	print('allele_ matrix:', allele_matrix)

	# cluster matrix should be empty
	assert len(cluster_matrix) == 0

	# allele matrix should be empty
	assert len(allele_matrix) == 0


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
	pruner = ReadSetPruning(readset, components, 2, 3, 3)
	cluster_matrix = pruner.get_cluster_matrix()
	allele_matrix = pruner.get_allele_matrix()
	print('cluster matrix:', cluster_matrix)
	print('allele matrix:', allele_matrix)

	expected_clusters = {
		0: [ ['Read 1', 'Read 2'], ['Read 3'] ],
		1: [ ['Read 4', 'Read 5'], ['Read 6'] ]
	}

	check_window_clustering(cluster_matrix, expected_clusters)
	consensus_clustering, optimal_partitioning = solve_MEC(cluster_matrix, 2)
	# TODO there are two unconnected blocks, they will be combined randomly
	expected_overall_clustering = sorted( [ ['Read 1', 'Read 2', 'Read 6'], ['Read 3', 'Read 4', 'Read 5'] ] )
	print('expected clusters: ', expected_overall_clustering)
	print('computed clusters: ', sorted(consensus_clustering))
	assert sorted(consensus_clustering) == expected_overall_clustering

# TODO: fails for number_of_clusters=2, since in the last window, it clusters the reads like this: {4,5}/{3}.
# 	Since similarity to read 3 is to different, the clusters are not combined.
#	However, solving MEC to get a consensus clustering, gives the correct solution
def test_clustering6():
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
	pruner = ReadSetPruning(readset, components, 2, 4, 3)
	cluster_matrix = pruner.get_cluster_matrix()
	allele_matrix = pruner.get_allele_matrix()
	print('cluster matrix: ', cluster_matrix)
	print('allele matrix:', allele_matrix)

	expected_clusters = {
		0:[ ['Read 1', 'Read 2'], ['Read 3', 'Read 4'] ],
		1:[ ['Read 2'], ['Read 3', 'Read 4', 'Read 5'] ]
	}

	check_window_clustering(cluster_matrix, expected_clusters)

	# get overall clustering solving MEC
	consensus_clustering, optimal_partitioning = solve_MEC(cluster_matrix, 2)
	expected_overall_clustering = sorted([ ['Read 1', 'Read 2'], ['Read 3', 'Read 4', 'Read 5'] ])
	print('expected clusters: ', expected_overall_clustering)
	print('computed clusters: ', sorted(consensus_clustering))
	assert sorted(consensus_clustering) == expected_overall_clustering


def test_clustering7():
	reads = """
                11
                1111
                  1111
                    11
		"""
	readset, positions, components = generate_input(reads)
	pruner = ReadSetPruning(readset, components, 2, 2, 2)
	cluster_matrix = pruner.get_cluster_matrix()
	allele_matrix = pruner.get_allele_matrix()
	print('cluster matrix: ', cluster_matrix)
	print('allele matrix:', allele_matrix)

	expected_clusters = {
		0:[ ['Read 1', 'Read 2'] ],
		1:[ ['Read 2', 'Read 3'] ],
		2:[ ['Read 3', 'Read 4'] ]
	}

	check_window_clustering(cluster_matrix, expected_clusters)

	# get overall clustering solving MEC
	consensus_clustering, optimal_partitioning = solve_MEC(cluster_matrix, 2)
	expected_overall_clustering = sorted([ [], ['Read 1', 'Read 2', 'Read 3', 'Read 4'] ])
	print('expected clusters: ', expected_overall_clustering)
	print('computed clusters: ', sorted(consensus_clustering))
	assert sorted(consensus_clustering) == expected_overall_clustering

