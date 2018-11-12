from tempfile import TemporaryDirectory
import os
import pysam
import math
import vcf

from pytest import raises
from whatshap.connect import run_connect

phasedblocks_onevariant = "tests/data/referencepanel/phasedblocks_onevariant.vcf"
reference_onevariant = "tests/data/referencepanel/reference_onevariant.vcf"


#tests one single variant
def test_one_variant():
	run_connect(variant_file=phasedblocks_onevariant, reference_file=reference_onevariant,
		output_file='/dev/null')
	
#tests variant file that contains the wrong chromosome (error expected)
def test_wrong_chromosome():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		with raises(SystemExit):
			run_connect(variant_file='tests/data/short-genome/wrongchromosome.vcf',	reference_file=reference_onevariant, output_file=outvcf)

#tests variant file using multiple chromosomes (error expected)
def test_multiple_chromosomes():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		with raises(SystemExit):
			run_connect(variant_file='tests/data/phased1.vcf',	reference_file=reference_onevariant, output_file=outvcf)	

#tests reference panel using multiple chromosomes (error expected)
def test_multiple_reference_chromosomes():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		with raises(SystemExit):
			run_connect(variant_file=phasedblocks_onevariant,	reference_file='tests/data/phased1.vcf', output_file=outvcf)


#tests variant file and reference panel that contain no common positions (error expected)
def test_no_variants():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		with raises(SystemExit):
			run_connect(variant_file='tests/data/referencepanel/phasedblocks_onevariant_no_common_positions.vcf',	reference_file=reference_onevariant, output_file=outvcf)

#tests writing the output file
def test_output():
	with TemporaryDirectory() as tempdir:
		outvcf = tempdir + '/output.vcf'
		run_connect(variant_file=phasedblocks_onevariant, reference_file=reference_onevariant,
		output_file=outvcf)
		lines = [l.strip('\n') for l in open(outvcf)]
		assert(len(lines)==7)