from nose.tools import raises
import math
from whatshap.core import ReadSet, PedigreeDPTable, Pedigree, NumericSampleIds, PhredGenotypeLikelihoods, GenotypeDPTable
from .phasingutils import string_to_readset, brute_force_phase
from .phasingutils import string_to_readset, string_to_readset_pedigree, brute_force_phase

def compare_to_expected(dp_forward_backward, positions, expected=None, genotypes=None):
	# check if computed likelihoods are equal to expected ones (if given)
	if not expected==None:
		for i in range(len(positions)):
			likelihoods = dp_forward_backward.get_genotype_likelihoods('individual0',i)
			print(likelihoods, expected[i], i)
			assert(likelihoods == expected[i])
		
	# check if likeliest genotype is equal to expected genotype
	for i in range(len(positions)):
		likelihoods = dp_forward_backward.get_genotype_likelihoods('individual0',i)
		# find likeliest genotype 
		max_val = -1
		max_index = -1
		for j in range(3):
			# likelihood should not be nan
			assert( not math.isnan(likelihoods[j]) )
			if likelihoods[j] > max_val:
				max_val = likelihoods[j]
				max_index = j
				
		print('genotype likelihoods for position',i, likelihoods, ' likeliest genotype: ',max_index)

		if not genotypes==None:
			assert(max_index==genotypes[i])


def test_genotyping_empty_readset():
	rs = ReadSet()
	genotypes = [1,1]
	recombcost = [1,1]
	numeric_sample_ids = NumericSampleIds()
	pedigree = Pedigree(numeric_sample_ids)
	genotype_likelihoods = [None, None]
	pedigree.add_individual('individual0', genotypes, genotype_likelihoods)
	dp_forward_backward = GenotypeDPTable(numeric_sample_ids,rs, recombcost, pedigree)
	

def check_genotyping_single_individual(reads, weights = None, expected = None, genotypes = None, scaling=None):
	# 0) set up read set
	readset = string_to_readset(reads, weights, scale_quality=scaling)
	positions = readset.get_positions()

	# 1) Genotype using forward backward algorithm
	recombcost = [1] * len(positions)
	numeric_sample_ids = NumericSampleIds()
	pedigree = Pedigree(numeric_sample_ids)
	genotype_likelihoods = [PhredGenotypeLikelihoods(0,0,0)] * len(positions)
	pedigree.add_individual('individual0', [1] * len(positions), genotype_likelihoods)
	dp_forward_backward = GenotypeDPTable(numeric_sample_ids, readset, recombcost,pedigree)

	# check the results
	compare_to_expected(dp_forward_backward, positions, expected, genotypes)
			
	# 2) Phase using PedMEC code for trios with two "empty" individuals (i.e. having no reads)
	recombcost = [1] * len(positions) # recombination costs 1, should not occur
	numeric_sample_ids = NumericSampleIds()
	pedigree = Pedigree(numeric_sample_ids)
	genotype_likelihoods = [PhredGenotypeLikelihoods(0,0,0)] * len(positions)
	pedigree.add_individual('individual0', [1] * len(positions), genotype_likelihoods)
	pedigree.add_individual('individual1', [1] * len(positions), genotype_likelihoods)
	pedigree.add_individual('individual2', [1] * len(positions), genotype_likelihoods)
	pedigree.add_relationship('individual0', 'individual1', 'individual2')
	dp_forward_backward = GenotypeDPTable(numeric_sample_ids,readset,recombcost,pedigree)
	
	# check the results
	compare_to_expected(dp_forward_backward, positions, expected, genotypes)
	
def test_geno_exact1() :
	reads = """
          11
           01
        """
	# as computed manually, with weight 10 for each position
	expected_likelihoods = [[0.05, 0.5, 0.45],[0.1323529411764706, 0.7352941176470589, 0.1323529411764706],[0.05, 0.5, 0.45]]
	genotypes = [1,1,1]
	check_genotyping_single_individual(reads,None,expected_likelihoods,genotypes,10)

def test_geno_exact2():
	reads = """
		11
		11
		"""
	weights = """
		11
		22
		"""
	# as computed manually, with given weights (x10)
	expected_likelihoods = [[0.0006656057777641766, 0.4062796462343544, 0.5930547479878814],[0.0006656057777641766, 0.4062796462343544, 0.5930547479878814]]
	genotypes = [2,2]
	check_genotyping_single_individual(reads,weights,expected_likelihoods,genotypes,10)

def test_geno_exact3():
	reads = """
          01
          11
        """
	# as computed manually, with weight 10 for each position
	expected_likelihoods = [[0.1493963782696177, 0.7012072434607646, 0.1493963782696177],[0.008551307847082495, 0.2987927565392354, 0.6926559356136821]]
	check_genotyping_single_individual(reads,None,expected_likelihoods,None,10)
   
def test_geno1():
	reads= """
	1111111111
	0000011111
	"""

	genotypes = [1,1,1,1,1,2,2,2,2,2]
	check_genotyping_single_individual(reads,None,None,genotypes,10)   
   
def test_geno2():
	reads = """
	101
	101
	101
	101
	100
	100
	100
	100
	"""
	
	genotypes = [2,0,1]
	check_genotyping_single_individual(reads,None,None,genotypes,10)

def test_geno3():
	reads = """
	011011
	110110
	110 10
	110110
	111110
	000 00
	01000 
	000010
	100100
	"""
	
	genotypes = [1,1,0,1,1,0]
	check_genotyping_single_individual(reads,None,None,genotypes,10)

##

def test_geno4():
	reads = """
	  1  11010
	  00 00101
	  001 01110
	   1    111
	"""

	check_genotyping_single_individual(reads,None,None,None,10)


def test_geno5():
	reads = """
	  0             0
	  110111111111
	  00100
	       0001000000
	       000
	        10100
	              101
	"""
	genotypes = [1,1,1,1,1,1,1,1,2,1,1,1,1,0,1]
	check_genotyping_single_individual(reads,None,None,genotypes,10)

def test_geno6():
	reads = """
		0100000000000
		0100010000000
		1110000000010
		0100000000000
		0101000001000
		0100010   000
		0 10000000100
		1111111011100
		0100111010011
		1111111000111
		1111110011111
		11110000  000
		1110011011111
		1111001011111
		0111111110  1
		"""
	genotypes = [1,2,1,1,1,1,1,0,1,1,1,1,1]
	check_genotyping_single_individual(reads,None,None,genotypes,60)

def test_geno7():
	reads = """
		111
		101
		111
		101
		010
		000
		010
		000
		"""
	genotypes = [1,1,1]
	check_genotyping_single_individual(reads,None,None,genotypes,60)

def test_geno8():
	reads = """
	11  
	11
	10
	"""
	genotypes = [2,1]
	check_genotyping_single_individual(reads, None,None,genotypes,10)

def test_geno9():
	reads = """
	111
	111
	010
	010
	   001
	   001
	   101
	   101
	"""
	genotypes = [1,2,1,1,0,2]
	check_genotyping_single_individual(reads,None,None,genotypes,10)

def test_weighted_genotyping1():
	reads = """
	  1  11010
	  00 00101
	  001 01110
	   1    111
	"""
	weights = """
	  2  13112
	  11 23359
	  223 56789
	   2    111
	"""
	
	genotypes = [1,1,1,1,1,1,2,1,1]
	check_genotyping_single_individual(reads, weights,None,genotypes,10)
	
def test_weighted_genotyping2():
	reads = """
	  111
	  101
	  111
	  101
	  010
	  000
	  010
	  000
	"""
	weights = """
	  999
	  999
	  999
	  999
	  999
	  999
	  999
	  999
	"""
	
	# second genotype should be twice as likely, since each (of the 4 different)
	# haplotype combinations has the same probability
	genotypes = [1,1,1]
	expected_likelihoods = [[0,1,0],[0.25,0.5,0.25],[0,1,0]]
	check_genotyping_single_individual(reads, weights,expected_likelihoods,genotypes,100)
	
def test_weighted_genotyping3():
	reads = """
		0 1
		 10
		 """
	weights = """
		999
		999
	"""
	expected_likelihoods = [[0.5,0.5,0],[0,0.5,0.5],[0,1,0]]
	check_genotyping_single_individual(reads,weights,expected_likelihoods, None, 500)		 
	
