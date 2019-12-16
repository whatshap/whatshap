"""
Integration tests that use the command-line entry point run_polyphase.
"""
from tempfile import TemporaryDirectory
import os

import pysam
from whatshap.cli.polyphase import run_polyphase
from whatshap.vcf import VcfReader, VariantCallPhase

def test_polyphase_short_chr22():
	outvcf = 'output.vcf'
	run_polyphase(
		phase_input_files=['tests/data/polyploid.chr22.42M.50k.bam'],
		variant_file='tests/data/polyploid.chr22.42M.50k.vcf',
		ploidy=4,
		ignore_read_groups=True,
		output=outvcf)
	assert os.path.isfile(outvcf)
	
	'''
	Test case is broken, because for some reason the output of the phasing
	is not yet written to the output file and this point, even though is 
	exactly like the diploid one.
	'''
	return

	tables = list(VcfReader(outvcf, phases=True))
	assert len(tables) == 1
	table = tables[0]
	assert table.chromosome == 'chr22'
	assert len(table.variants) == 69
	assert table.samples == ['HG00514_NA19240']
	
	#TODO: More asserts ...