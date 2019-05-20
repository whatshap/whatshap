from whatshap.core import ReadSet, compute_polyploid_genotypes
from whatshap.vcf import VcfGenotype
from whatshap.testhelpers import string_to_readset

def validate(reads, expected_genotypes, ploidy):
	readset = string_to_readset(reads, n_alleles=2)
	computed_genotypes = [VcfGenotype(gt) for gt in compute_polyploid_genotypes(readset, ploidy)]
	print('computed genotypes:', computed_genotypes)
	print('expected genotypes:', expected_genotypes)
	assert computed_genotypes == expected_genotypes

def test_simple1() :
	reads = """
	 1110
	 1100
	 1000
	"""
	expected_genotypes = [VcfGenotype([1,1,1]), VcfGenotype([1,1,0]), VcfGenotype([1,0,0]), VcfGenotype([0,0,0])];
	validate(reads, expected_genotypes, 3)

def test_simple2() :
	reads = """
	 1111
	 1011
	 1100
	 0100
	 0110
	 1111
	 1100
	"""
	expected_genotypes = [VcfGenotype([1,1,0]), VcfGenotype([1,1,1]), VcfGenotype([0,1,1]), VcfGenotype([0,0,1])];
	validate(reads, expected_genotypes, 3)

def test_simple3() :
	reads = """
	11110
	10000
	11010
	10010
	11010
	00111
	11111
	01110
	11110
	"""
	expected_genotypes = [VcfGenotype([1,1,1,0]), VcfGenotype([1,1,1,0]), VcfGenotype([1,1,0,0]), VcfGenotype([1,1,1,1]), VcfGenotype([0,0,0,1])]
	validate(reads, expected_genotypes, 4)
