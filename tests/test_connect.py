from tempfile import TemporaryDirectory
import os
import pysam
import math
import vcf

from pytest import raises
from whatshap.connect import included_variants
from whatshap.blockparsing import create_blocks, compute_referencepanel, compute_haplotypes, update_haplotypes
from whatshap.core import scoring_computation

phasedblocks_singletons = "tests/data/referencepanel/phasedblocks_singletons.vcf"
phasedblocks_normal = "tests/data/referencepanel/phasedblocks_normal.vcf"
reference_normal = "tests/data/referencepanel/reference_normal.vcf"
referencepanel_normal = "tests/data/referencepanel/reference_normal.txt"

phasedblocks_nested = "tests/data/referencepanel/phasedblocks_nested.vcf"
reference_nested = "tests/data/referencepanel/reference_nested.vcf"
referencepanel_nested = "tests/data/referencepanel/reference_nested.txt"

phasedblocks_contain_unphased = "tests/data/referencepanel/phasedblocks_contain_unphased.vcf"
referencepanel_contain_unphased = "tests/data/referencepanel/reference_contain_unphased.txt"

"""
tests if the block boundaries are created correctly according to the PS tag
"""
def test_creating_blocks_normal():
	variant_set = included_variants(reference_normal, phasedblocks_normal, False)
	assert(variant_set == {1,2,3,5,6,7,8,9,10})
	(blockends, haplo, intE) = create_blocks(phasedblocks_normal, "/dev/null", variant_set)
	assert (blockends == [3,6,8])
	assert(intE == [])
	assert(haplo == {1:('1', '0'),2:('1', '0'),3:('0', '1'), 5:('1', '0'), 6:('1','0'), 7:('0','1'),8:('0','1'),9:('1', '0'),10:('0','1') })

"""
tests if the reference panel is computed correctly 
"""
def test_computing_panel_normal():
	variant_set = included_variants(reference_normal, phasedblocks_normal, False)
	assert(variant_set == {1,2,3,5,6,7,8,9,10})
	
	compute_referencepanel(reference_normal, phasedblocks_normal, variant_set, referencepanel_normal)
	lines = [l.strip('\n') for l in open(referencepanel_normal)]
	#reference is supposed to contain 4 rows (haplotypes) and 9 columns (variants)
	assert(len(lines) == 4)
	assert(len(lines[0]) == 9)
	expected_lines = ["110101101","000011111","001010001","100000110"]
	for i in range(0,len(lines)):
		assert(lines[i] == expected_lines[i])

"""
input VCF: 3 blocks, no unphased variants
tests if the computation of the DP matrix and the backtracing find the correct best path pair and checks
the updating of the paths and corresponding haplotypes
"""
def test_path_computation():
	with TemporaryDirectory() as tempdir:
		haplotypefile = tempdir + "/haplotypes.txt"
		boundaryfile = tempdir + "/blockends.txt"
		costfile = tempdir + "/costs.txt"
		pathfile = tempdir + "/paths.txt"
		variant_set = included_variants(reference_normal, phasedblocks_normal, False)
		(blockends, haplo, internal_blockends) = create_blocks(phasedblocks_normal, boundaryfile, variant_set)
		assert(blockends == [3,6,8])
		assert(internal_blockends == [])
		compute_haplotypes(phasedblocks_normal, haplotypefile, variant_set)
		scoring_computation(haplotypefile, boundaryfile, referencepanel_normal, 1, costfile, pathfile)
		expected_haplos = ["110110010","001001101"]
		haplos = [l.strip('\n') for l in open(haplotypefile)]
		assert(len(haplos) == 2)
		for i in range(0,len(haplos)):
			assert(haplos[i] == expected_haplos[i])
		paths = [l.strip('\n') for l in open(pathfile)]
		assert(len(paths) == 2)
		expected_paths = ["0 0 0 0 2 2 2 3 3 ", "2 2 2 2 0 0 0 0 0 "]
		for i in range(0,len(paths)):
			assert(paths[i] == expected_paths[i])
		(newhaplo1, newhaplo2) = update_haplotypes(haplotypefile, blockends, pathfile, internal_blockends)
		expected_haplo1 = "110101101"
		expected_haplo2 = "001010010"
		assert(newhaplo1 == expected_haplo1 or newhaplo1 == expected_haplo2)
		assert(newhaplo2 == expected_haplo2 or newhaplo2 == expected_haplo1)

"""
pattern of the nested block IDs: 1 2 3 1 3 4
tests if all internal block boundaries are found correctly within the nested blocks, computes the DP matrix and
checks the correct path updating
"""
def test_update_paths_nestedblocks():
	variant_set = included_variants(reference_nested, phasedblocks_nested, False)
	assert(variant_set == {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,23,24,25,26,27,28,29,30,31,32})
	with TemporaryDirectory() as tempdir:
		boundaryfile = tempdir + "/blockends.txt"
		haplotypefile = tempdir + "/haplotypes.txt"
		costfile = tempdir + "/costs.txt"
		pathfile = "tests/data/referencepanel/paths.txt"
		(blockends, haplo, internal_blockends) = create_blocks(phasedblocks_nested, boundaryfile, variant_set)
		assert (blockends == [25,30])
		assert(internal_blockends == [[(0,5),(16,20)],[(6,10)]])
		compute_referencepanel(reference_nested, phasedblocks_nested, variant_set, referencepanel_nested)
		compute_haplotypes(phasedblocks_nested, haplotypefile, variant_set)
		expected_haplo1 = "1110111000101010110001111001001"		
		expected_haplo2 = "0001000111010101001110000110110"
		haplos = [l.strip('\n') for l in open(haplotypefile)]
		assert(len(haplos) == 2)
		assert(haplos[0] == expected_haplo1)	
		assert(haplos[1] == expected_haplo2)	
		scoring_computation(haplotypefile, boundaryfile, referencepanel_nested, 1, costfile, pathfile)
		(newhaplo1, newhaplo2) = update_haplotypes(haplotypefile, blockends, pathfile, internal_blockends)
		expected_new1 = "0001001000110101001110000101001"
		expected_new2 = "1110110111001010110001111010110"
		assert(newhaplo1 == expected_new1 or newhaplo1 == expected_new2)
		assert(newhaplo2 == expected_new2 or newhaplo2 == expected_new1)
"""
pattern of the nested block IDs: 1 2 3 2 4 1 5
tests if all internal block boundaries are found correctly within the nested blocks and
checks the correct path updating
"""
def test_update_paths_nestedblocks2():
	with TemporaryDirectory() as tempdir:
		haplotypefile = tempdir + "/haplotypes.txt"
		with open(haplotypefile,'w') as haplos:
			haplos.write("0011101111111110000101010000110"+'\n')
			haplos.write("1100010000000001111010101111001"+'\n')
		blockends = [25,30]
		internal_blockends = [[(4,7),(13,17)],[(8,12)],[(18,22)]]
		pathfile = tempdir + "/paths.txt"
		with open(pathfile,'w') as paths:	
			paths.write("1 1 1 1 2 2 2 2 2 2 2 2 2 3 3 3 3 3 1 1 1 1 1 4 4 4 2 2 2 2 2"+'\n')
			paths.write("2 2 2 2 3 3 3 3 3 3 3 3 3 4 4 4 4 4 2 2 2 2 2 2 2 2 4 4 4 4 4"+'\n')
		(newhaplo1, newhaplo2) = update_haplotypes(haplotypefile, blockends, pathfile, internal_blockends)
		expectedhaplo1 = "1100101100000110001010101100110"
		expectedhaplo2 = "0011010011111001110101010011001"
		assert(newhaplo1 == expectedhaplo1 or newhaplo1 == expectedhaplo2)
		assert(newhaplo2 == expectedhaplo2 or newhaplo2 == expectedhaplo1)
"""
input VCF contains only singletons (block size = 1)
tests if all positions are correctly detected as block boundaries
and checks the correct path computation and updating
"""
def test_singletons():
	with TemporaryDirectory() as tempdir:
		haplotypefile = tempdir + "/haplotypes.txt"
		boundaryfile = tempdir + "/blockends.txt"
		costfile = tempdir + "/costs.txt"
		pathfile = tempdir + "/paths.txt"
		variant_set = included_variants(reference_normal, phasedblocks_singletons, False)
		assert(variant_set == {1,2,3,4,5,6,7,8,9})
		compute_haplotypes(phasedblocks_singletons, haplotypefile, variant_set)
		expected_haplos = ["110110010","001001101"]
		haplos = [l.strip('\n') for l in open(haplotypefile)]
		assert(len(haplos) == 2)
		for i in range(0,len(haplos)):
			assert(haplos[i] == expected_haplos[i])
		(blockends, haplo, internal_blockends) = create_blocks(phasedblocks_singletons, boundaryfile, variant_set)
		assert(blockends == [0,1,2,3,4,5,6,7,8])
		assert(internal_blockends == [])
		scoring_computation(haplotypefile, boundaryfile, referencepanel_normal, 1, costfile, pathfile)		
		paths = [l.strip('\n') for l in open(pathfile)]
		assert(len(paths) == 2)
		expected_paths = ["0 0 0 0 2 2 2 3 3 ", "2 2 2 2 0 0 0 0 0 "]
		for i in range(0,len(paths)):
			assert(paths[i] == expected_paths[i])
		(newhaplo1, newhaplo2) = update_haplotypes(haplotypefile, blockends, pathfile, internal_blockends)
		expected_haplo1 = "001010010"
		expected_haplo2 = "110101101"
		assert(newhaplo1 == expected_haplo1 or newhaplo1 == expected_haplo2)
		assert(newhaplo2 == expected_haplo2 or newhaplo2 == expected_haplo1)

"""
input VCF contains unphased variant at position 3 (0-based)
tests if
	- the variant is correctly included within the panel and the haplotypes 
	- position 3 is detected as an internal block boundary
	- paths and haplotypes are correctly updated
"""
def test_contains_unphased():
	with TemporaryDirectory() as tempdir:
		boundaryfile = tempdir + "/blockends.txt"
		haplotypefile = tempdir + "/haplotypes.txt"
		costfile = tempdir + "/costs.txt"
		pathfile = tempdir + "/paths.txt"
		variant_set = included_variants(reference_normal, phasedblocks_contain_unphased, True)
		assert(variant_set == {1,2,3,4,5,6,7,8,9,10})
		compute_referencepanel(reference_normal, phasedblocks_contain_unphased, variant_set, referencepanel_contain_unphased)
		lines = [l.strip('\n') for l in open(referencepanel_contain_unphased)]
		assert(len(lines) == 4)
		assert(len(lines[0]) == 10)
		expected_lines = ["1100101101","0001011111","0010010001","1000000110"]
		for i in range(0,len(lines)):
			assert(lines[i] == expected_lines[i])
		(blockends, haplo, internal_blockends) = create_blocks(phasedblocks_contain_unphased, boundaryfile, variant_set)
		assert(blockends == [4,7,9])
		assert(internal_blockends == [[(3,3)]])
		compute_haplotypes(phasedblocks_contain_unphased, haplotypefile, variant_set)
		expected_haplos = ["1100110010","0011001101"]
		haplos = [l.strip('\n') for l in open(haplotypefile)]
		assert(len(haplos) == 2)
		for i in range(0,len(haplos)):
			assert(haplos[i] == expected_haplos[i])		
		scoring_computation(haplotypefile, boundaryfile, referencepanel_contain_unphased, 1, costfile, pathfile)		
		paths = [l.strip('\n') for l in open(pathfile)]
		assert(len(paths) == 2)
		expected_paths = ["0 0 0 0 0 2 2 2 3 3 ", "2 2 2 2 2 0 0 0 0 0 "]
		for i in range(0,len(paths)):
			assert(paths[i] == expected_paths[i])
		(newhaplo1, newhaplo2) = update_haplotypes(haplotypefile, blockends, pathfile, internal_blockends)
		expected_haplo1 = "0011010010"
		expected_haplo2 = "1100101101"
		assert(newhaplo1 == expected_haplo1 or newhaplo1 == expected_haplo2)
		assert(newhaplo2 == expected_haplo2 or newhaplo2 == expected_haplo1)

		
		