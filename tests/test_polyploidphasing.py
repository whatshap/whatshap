from whatshap.core import ReadSet, PedigreeDPTable, Pedigree, NumericSampleIds, PhredGenotypeLikelihoods
from whatshap.testhelpers import string_to_readset, brute_force_phase
from whatshap.phase import find_components
from whatshap.readsetpruning import ReadSetPruning
from whatshap.phase import find_components

def generate_cluster_input(reads, weights=None):
	readset = string_to_readset(reads, weights)
	positions = readset.get_positions()
	components = find_components(positions, readset)
	return readset, positions, components

def test_phase_empty_readset():
	rs = ReadSet()
	recombcost = [1,1]
	genotypes = [1,1]
	pedigree = Pedigree(NumericSampleIds(), 2)
	genotype_likelihoods = [None, None]
	pedigree.add_individual('individual0', genotypes, genotype_likelihoods)
	dp_table = PedigreeDPTable(rs, recombcost, pedigree, 2)
	superreads = dp_table.get_super_reads()

def rename_partitions(partitioning):
	last_index = 0
	mapping = {}
	for i in partitioning:
		if not i in mapping:
			mapping[i] = last_index
			last_index += 1
	return [mapping[i] for i in partitioning]

def compare_phasing_brute_force(superreads, cost, partition, readset, given_genotypes, ploidy, weights = None):
	"""Compares DPTable based phasing to brute force phasing and returns string representation of superreads."""
	assert len(superreads) == ploidy
	for i in range(1,ploidy):
		assert len(superreads[0]) == len(superreads[i])
	haplotypes = sorted(''.join(str(v.allele[0]) for v in sr) for sr in superreads)
	print(haplotypes, cost, partition)
	expected_cost, expected_partition, solution_count, expected_haplotypes = brute_force_phase(readset, ploidy, given_genotypes)
	print(haplotypes, expected_haplotypes)
	for sr in superreads:
		print(sr)
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

def check_phasing_single_individual(reads, genotypes, ploidy, weights = None, precomputed_readset = None):
	# 0) set up read set
	if precomputed_readset is None:
		readset = string_to_readset(reads, weights)
	else:
		readset = precomputed_readset
	positions = readset.get_positions()
	print(positions)
	# 1) Phase using PedMEC code for single individual
	for given_genotypes in [(True, genotypes), (False, None)]:
		recombcost = [1] * len(positions) # recombination costs 1, should not occur
		pedigree = Pedigree(NumericSampleIds(), ploidy)
		genotype_likelihoods = [None if given_genotypes[0] else PhredGenotypeLikelihoods([0] * (ploidy+1))] * len(positions)
		print(len(genotype_likelihoods), len(positions))
		pedigree.add_individual('individual0', genotypes, genotype_likelihoods)
		print("before DP table")
		dp_table = PedigreeDPTable(readset, recombcost, pedigree, ploidy, distrust_genotypes=not given_genotypes[0])
		print("after DP table")
		superreads, transmission_vector = dp_table.get_super_reads()
		print("after get_superreads")
		cost = dp_table.get_optimal_cost()
		# TODO: transmission vectors not returned properly, see issue 73
		assert len(set(transmission_vector)) == 1
		partition = dp_table.get_optimal_partitioning()
		compare_phasing_brute_force(superreads[0], cost, partition, readset, given_genotypes[1], ploidy, weights)

	# 2) Phase using PedMEC code for trios with two "empty" individuals (i.e. having no reads)
	for given_genotypes in [ (False, None), (True, genotypes)]:
		recombcost = [1] * len(positions) # recombination costs 1, should not occur
		pedigree = Pedigree(NumericSampleIds(), ploidy)
		genotype_likelihoods = [None if given_genotypes[0] else PhredGenotypeLikelihoods([0] * (ploidy+1))] * len(positions)
		pedigree.add_individual('individual0', genotypes, genotype_likelihoods)
		pedigree.add_individual('individual1', genotypes, genotype_likelihoods)
		pedigree.add_individual('individual2', genotypes, genotype_likelihoods)
		pedigree.add_relationship('individual0', 'individual1', 'individual2')
		dp_table = PedigreeDPTable(readset, recombcost, pedigree, ploidy, distrust_genotypes=not given_genotypes[0])
		cost = dp_table.get_optimal_cost()
		superreads, transmission_vector = dp_table.get_super_reads()
		assert len(set(transmission_vector)) == 1
		partition = dp_table.get_optimal_partitioning()
		compare_phasing_brute_force(superreads[0], cost, partition, readset, given_genotypes[1], ploidy, weights)


def test_phase_trivial() :
	reads = """
          111
          010
          101
        """
	check_phasing_single_individual(reads,[2,2,2], 3)

def test_phase1():
	reads = """
         11111
         00111
         11100
         11100
	"""
	check_phasing_single_individual(reads,[2,2,3,2,2], 3)

def test_phase2():
	reads = """
         1111111
         0101011
         0011110
         0111011
         1111101
	"""
	check_phasing_single_individual(reads, [1,2,3,3,2,3,2], 3)

# the long runtime of this tests is mainly due to the brute_force_phase method
def test_phase3():
	reads="""
        1111011
        0101110
        1110111
        1010111
        1010111
        0000001
        0000000
	"""
	check_phasing_single_individual(reads, [2,2,2,2,2,3,3], 4)

def test_phase4():
	reads="""
        00000
        00010
        00101
        01101
        10001
        01110
	"""
	check_phasing_single_individual(reads, [1,1,2,1,1], 4)


# this takes 15 min to be solved ...
#def test_phase4():
#	reads="""
#        11110000
#        10101010
#        01010101
#        00001111
#        11111111
#        00000000
#        11001100
#        00110011 
#	"""
#	check_phasing_single_individual(reads, [4,4,4,4,3,3,3,3], 8)
#
#def test_phase5():
#	reads="""
#        11111111
#        11111111
#        01010101
#        01100101
#        10101010
#        10101011
#        00000000
#        00001111
#	"""
#	check_phasing_single_individual(reads, [2,2,2,2,3,3,3,4], 7)

def test_clustered_phasing1():
	reads = """
                111
                111
                010
                000
                100
                100
		"""
	weights = """
                111
                111
                111
                222
                333
                111
		"""
	# first cluster some reads
	readset, positions, components = generate_cluster_input(reads, weights)
	pruner = ReadSetPruning(readset, components, 3, False)
	pruned_readset = pruner.get_pruned_readset()
	print("pruned readset: ", pruned_readset)
	check_phasing_single_individual(None, [1,1,1], 3,  precomputed_readset = pruned_readset)

def test_clustered_phasing2():
	reads = """
                11111111
                01111101
                00011000
                00001000
                10101010
                10101011
                11100001
                11100000
		"""
	weights = """
                23543243
                11123222
                22211111
                24333231
                12121222
                11111222
                12122212
                11133332
		"""
	readset, positions, components = generate_cluster_input(reads,weights)
	pruner = ReadSetPruning(readset, components, 4, False)
	pruned_readset = pruner.get_pruned_readset()
	print("pruned readset: ", pruned_readset)
	check_phasing_single_individual(None, [3,2,3,1,3,1,2,1], 4, precomputed_readset = pruned_readset)

def test_clustered_phasing3():
	reads = """
                111111111111111111111
                11100101000100001001
                00011000000011000000
                  01100000001100110000
                      1111111111101111111
                      01000100011001101
		"""
	readset, positions, components = generate_cluster_input(reads)
	pruner = ReadSetPruning(readset, components, 3, False)
	pruned_readset = pruner.get_pruned_readset()
	print('pruned readset: ', pruned_readset)
	check_phasing_single_individual(None, [2,2,2,2,2,2,1,2,1,1,1,2,2,2,1,1,2,1,1,2,2,1,2,1,1], 3, precomputed_readset = pruned_readset)
