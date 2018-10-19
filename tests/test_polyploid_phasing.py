from whatshap.core import ReadSet, PedigreeDPTable, Pedigree, NumericSampleIds, GenotypeLikelihoods, Genotype
from whatshap.testhelpers import string_to_readset, brute_force_phase
from whatshap.phase import find_components
from whatshap.clusterediting import print_readset
from whatshap.readsetpruning import ReadSetPruning
from whatshap.matrixtransformation import MatrixTransformation

def create_genotype_vector(ploidy, expected_genotypes):
	genotype_vector = []
	for gt in expected_genotypes:
		n_ref_alleles = ploidy - gt
		n_alt_alleles = gt
		alleles = [0]*n_ref_alleles + [1]*n_alt_alleles
		genotype_vector.append(Genotype(alleles))
	return genotype_vector

def solve_MEC(cluster_matrix, ploidy, cluster_counts):
	"""
	finds the best consensus clustering solving MEC.
	"""
	print_readset(cluster_matrix)
	numeric_sample_ids = NumericSampleIds()
	pedigree = Pedigree(numeric_sample_ids, ploidy)
	windows = cluster_matrix.get_positions()
	# TODO compute number of genotypes
	genotypes = [Genotype([i for i in range(0,ploidy)])] * len(windows)
	print(genotypes, cluster_counts)
	pedigree.add_individual('0',genotypes,[GenotypeLikelihoods(ploidy, ploidy, [])]*len(windows))
	dp_table = PedigreeDPTable(cluster_matrix, [1]*len(windows), pedigree, ploidy, distrust_genotypes=False, allele_counts=cluster_counts)
	result = []
	print('cost: ', dp_table.get_optimal_cost())
	for i in range(ploidy):
		result.append([])
	optimal_partitioning = dp_table.get_optimal_partitioning()
	print(optimal_partitioning)
	for i,partition in enumerate(optimal_partitioning):
		result[partition].append(cluster_matrix[i].name)
	for i in range(len(result)):
		result[i] = sorted(result[i])
	return sorted(result), optimal_partitioning

def derive_haplotypes(reads, positions, ploidy, given_genotypes, precomputed_partitioning):
	"""
	given reads and a partitioning, computes the best allele assignment.
	"""
	numeric_sample_ids = NumericSampleIds()
	print_readset(reads)
	pedigree = Pedigree(numeric_sample_ids, ploidy)
	genotype_likelihoods = [None if given_genotypes[0] else GenotypeLikelihoods(ploidy, 2, [0] * (ploidy+1))] * len(positions)
	genotypes = given_genotypes[1] if given_genotypes[0] else [Genotype([])]*len(positions)
	pedigree.add_individual('0', genotypes, genotype_likelihoods)
	dp_table = PedigreeDPTable(reads, [1]*len(positions), pedigree, ploidy, distrust_genotypes=not given_genotypes[0], positions=positions, precomputed_partitioning=precomputed_partitioning)
	print('cost 2: ', dp_table.get_optimal_cost())
	return dp_table.get_super_reads()[0], dp_table.get_optimal_cost()

def reorder_optimal_partitioning(readset1, partitioning, readset2):
	"""
	order the partitioning according to the read order in readset2.
	"""
	# map read id to partition (since sorting of reads can be different)
	read_to_partition = {}
	for i,read in enumerate(readset1):
		read_to_partition[read.name] = partitioning[i]
	return [ read_to_partition[r.name] for r in readset2]

def rename_partitions(partitioning):
	last_index = 0
	mapping = {}
	for i in partitioning:
		if not i in mapping:
			mapping[i] = last_index
			last_index += 1
	return [mapping[i] for i in partitioning]

def compare_phasing_brute_force(superreads, cost, partition, readset, given_genotypes, ploidy, weights = None):
	"""
	Compares DPTable based phasing to brute force phasing and returns string representation of superreads.
	"""
	assert len(superreads) == ploidy
	for i in range(1,ploidy):
		assert len(superreads[0]) == len(superreads[i])
	haplotypes = sorted(''.join(str(v.allele) for v in sr) for sr in superreads)
	expected_cost, expected_partition, solution_count, expected_haplotypes = brute_force_phase(readset, ploidy, 2, given_genotypes)
	print('Partition:', partition, rename_partitions(partition))
	print('Expected: ', expected_partition, rename_partitions(expected_partition))
	print('Haplotypes:')
	for h in haplotypes:
		print(h)
	print('Expected Haplotypes:')
	for h in expected_haplotypes:
		print(h)
	print('Cost:', cost)
	print('Expected cost:', expected_cost)
	assert (rename_partitions(partition) == rename_partitions(expected_partition))
	assert solution_count == 1
	assert cost == expected_cost
	assert(sorted(haplotypes) == sorted(expected_haplotypes))


def check_phasing_single_individual(reads, genotypes, ploidy, reads_per_window, variants_per_window, algorithm, weights = None):
	# 0) set up read set
	readset = string_to_readset(reads, weights)
	positions = readset.get_positions()
	components = find_components(positions, readset)

	# 1) construct partitioning of the reads
	if algorithm == 'windowphase':
		pruner = ReadSetPruning(readset, components, ploidy, reads_per_window, variants_per_window)
		allele_matrix = pruner.get_allele_matrix()
		cluster_matrix = pruner.get_cluster_matrix()
	elif algorithm == 'clusterediting':
		pruner = MatrixTransformation(readset, components, ploidy, 0.1, variants_per_window)
		allele_matrix = readset
		cluster_matrix = pruner.get_transformed_matrix()
	else:
		assert(False)
	cluster_counts = pruner.get_cluster_counts()
	print('cluster matrix:', cluster_matrix)
	print('allele_counts:', cluster_counts)
	print('allele matrix:', allele_matrix)
	# combine columnwise partitionings
	consensus_clustering, optimal_partitioning = solve_MEC(cluster_matrix, ploidy, cluster_counts)
	print('optimal_partitioning:', optimal_partitioning)

	# 2) given the partitioning of reads, derive the optimal haplotypes
	ordered_partitioning = reorder_optimal_partitioning(cluster_matrix, optimal_partitioning, allele_matrix)
	print('ordered_partitioning:', ordered_partitioning)
	for given_genotypes in [ (True, genotypes), (False, None)]:
		superreads, optimal_cost = derive_haplotypes(allele_matrix, positions, ploidy, given_genotypes, ordered_partitioning)
		final_partitioning = reorder_optimal_partitioning(allele_matrix, ordered_partitioning, readset)
		compare_phasing_brute_force(superreads[0], optimal_cost, final_partitioning, readset, given_genotypes[1], ploidy, weights)

def test_diploid_phase1():
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
	for algorithm in ['windowphase', 'clusterediting']:
		genotypes = create_genotype_vector(2, [1,1,2,1,2,2,1,1])
		check_phasing_single_individual(reads, genotypes, 2, 3, 6, algorithm)

def test_diploid_phase2():
	reads = """
                11111
                111110
                110111
                 010010
                  10010
                    0100
                    01001000
                     1110111
                        1111
		"""
	for algorithm in ['windowphase', 'clusterediting']:
		genotypes = create_genotype_vector(2, [1,1,2,1,1,2,1,1,1,1,1,1])
		check_phasing_single_individual(reads, genotypes, 2, 3, 3, algorithm)

def test_diploid_phase3():
	reads = """
	  11 11010
	  000001010
	  001001110
	   11   111
	"""
	for algorithm in ['windowphase', 'clusterediting']:
		genotypes = create_genotype_vector(2, [1,1,1,1,1,1,1,1,1])
		check_phasing_single_individual(reads, genotypes, 2, 3, 3, algorithm)

def test_diploid_phase4():
	reads = """
         1111111
         0000000
            111111111
              0000000
	"""
	for algorithm in ['windowphase', 'clusterediting']:
		genotypes = create_genotype_vector(2, [1,1,1,1,1,1,1,1,1,1,1,1])
		check_phasing_single_individual(reads, genotypes, 2, 2, 3, algorithm)

# TODO approach to contructing partitioning works badly in such examples since intersection
# of variants covered by all reads in a window is small
# ==> several unconnected blocks and many unphased variants result
#def test_diploid_phase5():
#	reads = """
#	  0             0
#	  110111111111
#	  00100
#	       0001000000
#	       000
#	        10100
#	              101
#	"""
#	genotypes = create_genotype_vector(2, [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
#	check_phasing_single_individual(reads, genotypes, 2, 4, 3, 'clusterediting')


# TODO: cluster editing always puts Read 1 and Read 3 in the same cluster since they are too similar
def test_polyploid_phase1():
	reads = """
          111
          010
          101
	"""
	for algorithm in ['windowphase']:
		genotypes = create_genotype_vector(3, [2,2,2])
		check_phasing_single_individual(reads, genotypes, 3, 3, 3, algorithm)

def test_polyploid_phase2():
	reads = """
         11111
         00111
         11100
         11100
	"""
	for algorithm in ['windowphase', 'clusterediting']:
		genotypes = create_genotype_vector(3, [2,2,3,2,2])
		check_phasing_single_individual(reads, genotypes, 3, 3, 5, algorithm)

def test_polyploid_phase3():
	reads = """
         1111111
         0101011
         0011110
         0111011
         1111101
	"""
	for algorithm in ['windowphase', 'clusterediting']:
		genotypes = create_genotype_vector(3, [1,2,3,3,2,3,2])
		check_phasing_single_individual(reads, genotypes, 3, 3, 5, algorithm)

# TODO: windowphase: connection reads 1 and 5,6 are never in the same window.
#	they will be put in the same set, since they are never compared
#	and the others do not give any 'hints' on whether to cluster
#	1 / 5,6 / ... or  1,5,6 / ..
#	leads to wrong partitioning and much higher MEC score
#	clusterediting: reads 1 and 3,4,5 are put in the same cluster,
#	since they are too similar
#def test_polyploid_phase4():
#	reads="""
#        1111011
#        0101110
#        1110111
#        1010111
#        1010111
#        0000001
#        0000000
#	"""
#	genotypes = create_genotype_vector(4, [2,2,2,2,2,3,3])
#	check_phasing_single_individual(reads, genotypes, 4, 5, 5, 'clusterediting')

def test_polyploid_phase5():
	reads="""
        00000
        00010
        00101
        01101
        10001
        01110
	"""
	for algorithm in ['windowphase', 'clusterediting']:
		genotypes =  create_genotype_vector(4, [1,1,2,1,1])
		check_phasing_single_individual(reads, genotypes, 4, 6, 5, algorithm)


