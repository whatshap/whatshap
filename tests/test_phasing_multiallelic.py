from whatshap.core import ReadSet, PedigreeDPTable, Pedigree, NumericSampleIds, GenotypeLikelihoods, Genotype
from whatshap.testhelpers import string_to_readset, brute_force_phase
import math

def binomial_coeff(x,y):
	if y == x:
    		return 1
	elif y == 1:
		return x
	elif y > x:
		return 0
	else:
		a = math.factorial(x)
		b = math.factorial(y)
		c = math.factorial(x-y)
		div = a // (b * c)
	return div


def rename_partitions(partitioning):
	last_index = 0
	mapping = {}
	for i in partitioning:
		if not i in mapping:
			mapping[i] = last_index
			last_index += 1
	return [mapping[i] for i in partitioning]


def compare_phasing_brute_force(superreads, cost, partition, readset, given_genotypes, ploidy, n_alleles, weights = None):
        """
        Compares DPTable based phasing to brute force phasing and returns string representation of superreads.
        """
        assert len(superreads) == ploidy
        for i in range(1,ploidy):
                assert len(superreads[0]) == len(superreads[i])
        haplotypes = sorted(''.join(str(v.allele) for v in sr) for sr in superreads)
        expected_cost, expected_partition, solution_count, expected_haplotypes = brute_force_phase(readset, ploidy, n_alleles, given_genotypes)
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


def check_phasing_single_individual(reads, ploidy, n_alleles, given_genotypes, weights = None):
	# 0) set up read set
	readset = string_to_readset(reads, n_alleles=n_alleles, w=weights)
	positions = readset.get_positions()
	n_genotypes = binomial_coeff(ploidy + n_alleles - 1, n_alleles - 1)

	# 1) Phase using PedMEC code for single individual
	for gt_given in [False, True]:
		recombcost = [1] * len(positions) # recombination costs 1, should not occur
		pedigree = Pedigree(NumericSampleIds(), ploidy)
		genotype_likelihoods = [None if gt_given else GenotypeLikelihoods(ploidy,n_alleles,[0]*n_genotypes)] * len(positions)
		genotypes = given_genotypes if gt_given else [Genotype([])] * len(positions)
		pedigree.add_individual('individual0', genotypes, genotype_likelihoods) # all genotypes heterozygous
		print('test_phasing: ', readset, pedigree)
		dp_table = PedigreeDPTable(readset, recombcost, pedigree, ploidy, distrust_genotypes=not gt_given, allele_counts = [n_alleles] * len(positions))
		superreads, transmission_vector = dp_table.get_super_reads()
		print("after get_superreads")
		cost = dp_table.get_optimal_cost()
		assert len(set(transmission_vector)) == 1
		partition = dp_table.get_optimal_partitioning()
		allowed_genotypes = given_genotypes if gt_given else None
		print("allowed genotypes: ", allowed_genotypes, "distrust genotypes: ", not gt_given)
		compare_phasing_brute_force(superreads[0], cost, partition, readset, allowed_genotypes, ploidy, n_alleles, weights)


#def test_phase1():
#	reads = """
#	 10
#	 010
#	 010
#	"""
#	check_phasing_single_individual(reads)

def test_multiallelic_phase1():
	reads = """
	 10
	 020
	 021
	"""
	genotypes = [Genotype([0,1]), Genotype([0,2]), Genotype([1,1])]
	check_phasing_single_individual(reads, 2, 4, genotypes)


def test_multiallelic_phase2():
	reads = """
	  1230
	  1230
	  012301
	  0123
	    301232
	     0 232
	     30223
	       223
	        23
	"""
	genotypes = [Genotype([0,1]),Genotype([1,2]),Genotype([2,3]),Genotype([0,3]),Genotype([0,1]),Genotype([2,2]),Genotype([2,3]),Genotype([2,3])]
	check_phasing_single_individual(reads, 2, 4, genotypes)

def test_multiallelic_phase3():
	reads = """
	  1230
	  1230
	  012301
	  0123
	    301232
	     0 232
	     30223
	       223
	        23
	"""
	genotypes = [Genotype([0,1]),Genotype([1,2]),Genotype([2,3]),Genotype([0,3]),Genotype([0,1]),Genotype([0,0]),Genotype([2,3]),Genotype([2,3])]
	check_phasing_single_individual(reads, 2, 4, genotypes)

def test_multiallelic_phase4():
	reads = """
	 011
	  0213
	  021
	  101
	   002
	    12
	"""

	weights = """
	 114
	  1311
	  132
	  132
	   321
	    11
        """
	genotypes = [Genotype([0,1]), Genotype([0,1]), Genotype([0,2]), Genotype([1,1]), Genotype([2,3])]
	check_phasing_single_individual(reads, 2, 4, genotypes, weights)

def test_multiallelic_phase5():
	reads =  """
	 00000
	 00010
	 11101
	 11111
	 22022
	 22212
	"""
	genotypes = [Genotype([0,1,2])] * 5
	check_phasing_single_individual(reads, 3, 3, genotypes)	
