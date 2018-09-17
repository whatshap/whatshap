from whatshap.core import ReadSet, PedigreeDPTable, Pedigree, NumericSampleIds, PhredGenotypeLikelihoods, Genotype
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
	print(position_to_clusters)

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
	# TODO number of possible genotypes
	# allowed genotypes (require different allele for each partition, e.g. gt 0/1/2/3 for ploidy=4)
	genotypes = [Genotype([i for i in range(0,ploidy)])] * len(windows)
	pedigree.add_individual('0', genotypes, [PhredGenotypeLikelihoods(ploidy, ploidy,[0]*(ploidy+1))]*len(windows))
	# n_alleles per window == ploidy
	allele_counts = [ploidy] * len(windows)
	dp_table = PedigreeDPTable(cluster_matrix, [1]*len(windows), pedigree, ploidy, distrust_genotypes=True, allele_counts=allele_counts)
	result = []
	for i in range(ploidy):
		result.append([])
	optimal_partitioning = dp_table.get_optimal_partitioning()
	for i,partition in enumerate(optimal_partitioning):
		result[partition].append(cluster_matrix[i].name)
	for i in range(len(result)):
		result[i] = sorted(result[i])
	return sorted(result), optimal_partitioning

#def best_allele_assignment(reads, positions, ploidy, genotypes, precomputed_partitioning):
#	"""
#	given reads and a partitioning, computes the best allele assignment.
#	"""
#	numeric_sample_ids = NumericSampleIds()
#	pedigree = Pedigree(numeric_sample_ids, ploidy)
#	pedigree.add_individual('0', genotypes,	[PhredGenotypeLikelihoods([0]*(ploidy+1))]*len(genotypes))
#	dp_table = PedigreeDPTable(reads, [1]*len(positions), pedigree, ploidy, False, positions, precomputed_partitioning)
#	print('SUPERREADS:', dp_table.get_super_reads()[0][0], dp_table.get_optimal_cost())


def reorder_optimal_partitioning(readset1, partitioning, readset2):
	"""
	order the partitioning according to the read order in readset2.
	"""
	# map read id to partition (since sorting of reads can be different)
	read_to_partition = {}
	for i,read in enumerate(readset1):
		read_to_partition[read.name] = partitioning[i]
	return [ read_to_partition[r.name] for r in readset2]

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

	# get overall clustering solving MEC	
	consensus_clustering, optimal_partitioning = solve_MEC(cluster_matrix, 2)
	assert sorted(consensus_clustering) == sorted([ ['Read 1', 'Read 2'], ['Read 3', 'Read 4', 'Read 5'] ])
#	genotypes = [1,1,2,1,2,2,1,1]
#	ordered_partitioning = reorder_optimal_partitioning(cluster_matrix, optimal_partitioning, readset)
#	best_allele_assignment(readset, positions, 2, genotypes, ordered_partitioning)

# TODO: disadvantage of the approach: Since only reads in a window are considered, and not all covering a region at the same time
#	the clustering misses some of the relationships between the reads. Here Read3 and Read8 are not combined, because they
#	are never within the same window (and no other reads in between relate to them)
# TODO: the implementation does not yet work for ploidy > 2
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
	cluster_matrix = pruner.get_cluster_matrix()
	allele_matrix = pruner.get_allele_matrix()
	print('cluster matrix', cluster_matrix)
	print('allele_matrix', allele_matrix)

	expected_clusters = {
		0: [ ['Read 1', 'Read 2'], ['Read 3'], ['Read 4'], ['Read 5'] ],
		1: [ ['Read 2'], ['Read 3'], ['Read 4'], ['Read 5', 'Read 6'] ],
		2: [ ['Read 3'], ['Read 4', 'Read 7'], ['Read 5', 'Read 6'] ],
		3: [ ['Read 4', 'Read 7'], ['Read 5', 'Read 6'], ['Read 8'] ]
	}

	check_window_clustering(cluster_matrix, expected_clusters)

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
	assert sorted(consensus_clustering) == sorted( [ ['Read 1', 'Read 2', 'Read 6'], ['Read 3', 'Read 4', 'Read 5'] ] )
#	genotypes = [1,1,2,1,1,1,2,1,1,0,1,1]
#	ordered_partitioning = reorder_optimal_partitioning(cluster_matrix, optimal_partitioning, readset)
#	best_allele_assignment(readset, positions, 2, genotypes, ordered_partitioning)

# TODO: adjust to new version, requires multiallelic MEC to work
# TODO: same issue as before: especially as ploidy increases, the number of reads per window needs to be large in order to not miss to many
#	relationships between reads.
#def test_clustering6():
#	reads = """
#                1010001001101010001010100010001111111010100100
#                1110100010001000101000000101000101100100011101000
#                000010010101101100010001111000011000011000010000011
#                000011111111100000100000011100001111111111011111001111
#                101000100000000010001111100000010000100011100011111010
#                101000000000000010001111101000000010000011100011111010
#                010101010000000011100111111100011100001111110001011110
#                010101010000000000000111111111111100001111110001111010
#                000011000111000111001001111100001010101010001001001001
#                001010100101010000000100101101101110111100011001011000
#                001000100101010001100000101111101110111100011111011000
#                   011110111100010000000000000001111111101111111001111
#                    00010101111100010001111000000010111000000010011111
#                    000010001010101000010101010100010101010100100111
#                      000101000111001001111000101010101010001001001
#                        1010101000100010100100111111101010010010001011
#		"""
#	readset, positions, components = generate_input(reads)
#	expected_clusters = [['Read 1', 'Read 16'],['Read 14', 'Read 2'],['Read 13', 'Read 3'], ['Read 12', 'Read 4'],['Read 5', 'Read 6'],['Read 7', 'Read 8'],
#				['Read 15', 'Read 9'], ['Read 10', 'Read 11'] ]
#	pruner = ReadSetPruning(readset, components, 8, 20, 5)
#	computed_clusters = pruner.get_clusters()
#	print("computed: ", computed_clusters)
#	print("expected: ", expected_clusters)
#	assert computed_clusters == expected_clusters


# TODO: fails for number_of_clusters=2, since in the last window, it clusters the reads like this: {4,5}/{3}.
# 	Since similarity to read 3 is to different, the clusters are not combined.
#	However, solving MEC to get a consensus clustering, gives the correct solution
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
	pruner = ReadSetPruning(readset, components, 2, 3, 3)
	cluster_matrix = pruner.get_cluster_matrix()
	allele_matrix = pruner.get_allele_matrix()
	print('cluster matrix: ', cluster_matrix)
	print('allele matrix:', allele_matrix)

	expected_clusters = {
		0:[ ['Read 1', 'Read 2'], ['Read 3'] ],
		1:[ ['Read 2'], ['Read 3', 'Read 4'] ],
		2:[ ['Read 3'], ['Read 4', 'Read 5'] ]
	}

	check_window_clustering(cluster_matrix, expected_clusters)

	# get overall clustering solving MEC
	consensus_clustering, optimal_partitioning = solve_MEC(cluster_matrix, 2)
	assert sorted(consensus_clustering) == sorted([ ['Read 1', 'Read 2'], ['Read 3', 'Read 4', 'Read 5'] ])
#	genotypes = [1,1,2,1,2,1]
#	ordered_partitioning = reorder_optimal_partitioning(cluster_matrix, optimal_partitioning, readset)
#	best_allele_assignment(readset, positions, 2, genotypes, ordered_partitioning)

def test_clustering8():
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
	assert sorted(consensus_clustering) == sorted([ [], ['Read 1', 'Read 2', 'Read 3', 'Read 4'] ])
#	genotypes = [2,2,2,2,2,2]
#	ordered_partitioning = reorder_optimal_partitioning(cluster_matrix, optimal_partitioning, readset)
#	best_allele_assignment(readset, positions, 2, genotypes, ordered_partitioning)
