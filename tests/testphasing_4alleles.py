from nose.tools import raises
from whatshap.core import ReadSet, PedigreeDPTable, Pedigree, NumericSampleIds, PhredGenotypeLikelihoods
from .phasingutils import string_to_readset, brute_force_phase


def phase_MAV(reads, n_alleles, all_het, genos, genotypes, weights =  None):
	readset = string_to_readset(reads, n_alleles)
	positions = readset.get_positions()
	for all_heterozygous in all_het:
		recombcost = [1] * len(positions) # recombination costs 1, should not occur
		pedigree = Pedigree(NumericSampleIds())
		genotype_likelihoods = [None if all_heterozygous else PhredGenotypeLikelihoods(genos)] * len(positions)
		pedigree.add_individual('individual0', genotypes, genotype_likelihoods) # all genotypes heterozygous
		dp_table = PedigreeDPTable(readset, recombcost, pedigree, distrust_genotypes=not all_heterozygous)
		superreads_list, transmission_vector = dp_table.get_super_reads()
		cost = dp_table.get_optimal_cost()
	return superreads_list, transmission_vector, cost
      
def assert_haplotypes(superreads_list, all_expected_haplotypes, length):
	for superreads, expected_haplotypes in zip(superreads_list, all_expected_haplotypes):
		assert len(superreads) == 2
		assert len(superreads[0]) == len(superreads[1]) == length
		haplotypes = tuple(sorted(''.join(str(v.allele) for v in sr) for sr in superreads))
		assert (haplotypes == (expected_haplotypes[0], expected_haplotypes[1])) or (haplotypes == (expected_haplotypes[1], expected_haplotypes[0]))

def test_phase1():
	reads = """
	 10
	 020
	 021
	"""
	superreads_list, transmission_vector, cost = phase_MAV(reads, 4, [True], [5,5,5,5,5,5,5,5,5,5], [1,2,2])
	assert cost == 1
	all_expected_haplotypes = [
		('02-2','10-2'),
	]
	#('020','102'),
	assert_haplotypes(superreads_list, all_expected_haplotypes, 3)

def test_phase2():
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
	superreads_list, transmission_vector, cost = phase_MAV(reads, 4, [True, False], [5,5,5,5,5,5,5,5,5,5], [1,3,5,3,1,4,5,5])
	assert cost == 41
	all_expected_haplotypes = [
		('01230223','12301232'),
	]
	assert_haplotypes(superreads_list, all_expected_haplotypes, 8)

def test_phase2_flip_gl():
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
	superreads_list, transmission_vector, cost = phase_MAV(reads, 4, [False], [5,5,5,5,0,5,5,5,5,5,5], [1,3,5,3,1,0,5,5])
	assert cost == 17
	all_expected_haplotypes = [
		('1-2331223','3-2113221'),
	]
	#('11331223','33113221'),
	assert_haplotypes(superreads_list, all_expected_haplotypes, 8)


def test_phase2_flip_gt():
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
	superreads_list, transmission_vector, cost = phase_MAV(reads, 4, [False], [5,5,5,5,5,5,5,5,5,5,5], [1,3,5,3,1,0,5,5])
	assert cost == 41
	all_expected_haplotypes = [
		('01230223','12301232'),
	]
	assert_haplotypes(superreads_list, all_expected_haplotypes, 8)
	
def test_phase2_flip_gt_allzeros():
	reads = """
	  1230
	  1230
	  012301
	  0123
	    301232
	     1 232
	     30223
	       223
	        23
	"""
	superreads_list, transmission_vector, cost = phase_MAV(reads, 4, [False], [0,0,0,0,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,0])
	assert cost == 2
	all_expected_haplotypes = [
		('01230223','12301232'),
	]
	assert_haplotypes(superreads_list, all_expected_haplotypes, 8)
