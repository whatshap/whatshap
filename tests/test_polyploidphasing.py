from whatshap.core import ReadSet, PedigreeDPTable, Pedigree, NumericSampleIds, PhredGenotypeLikelihoods
from whatshap.testhelpers import string_to_readset, brute_force_phase


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
	haplotypes = sorted(''.join(str(v.allele) for v in sr) for sr in superreads)
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

def check_phasing_single_individual(reads, genotypes, ploidy, weights = None):
	# 0) set up read set
	readset = string_to_readset(reads, weights)
	positions = readset.get_positions()

	# 1) Phase using PedMEC code for single individual
	for given_genotypes in [ (False, None), (True, genotypes)]:
		recombcost = [1] * len(positions) # recombination costs 1, should not occur
		pedigree = Pedigree(NumericSampleIds(), ploidy)
		genotype_likelihoods = [None if given_genotypes[0] else PhredGenotypeLikelihoods([0] * (ploidy+1))] * len(positions)
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
