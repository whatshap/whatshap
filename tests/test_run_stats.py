"""
Tests for 'whatshap stats'
"""

from tempfile import TemporaryDirectory
from collections import namedtuple
from whatshap.stats import run_stats, run_stats_bed

def test_stats1():
	with TemporaryDirectory() as tempdir:
		outtsv = tempdir + '/output.tsv'
		run_stats(ploidy=2, vcf='tests/data/phased1.vcf', tsv=outtsv, sample='sample1', chr_lengths='tests/data/chr-lengths.txt')
		lines = [l.split('\t') for l in open(outtsv)]
		assert len(lines) == 4
		Fields = namedtuple('Fields', [f.strip('#\n') for f in lines[0]])
		entry_chrA, entry_chrB, entry_all = [Fields(*l) for l in lines[1:]]

		assert entry_chrA.chromosome == 'chrA'
		assert entry_chrA.variants == '8'
		assert entry_chrA.phased == '7'
		assert entry_chrA.unphased == '1'
		assert entry_chrA.blocks == '2'
		assert entry_chrA.variant_per_block_sum == '7'
		assert entry_chrA.bp_per_block_sum == '551'
		assert entry_chrA.block_n50[:-1] == '101'

		assert entry_chrB.chromosome == 'chrB'
		assert entry_chrB.variants == '2'
		assert entry_chrB.phased == '2'
		assert entry_chrB.unphased == '0'
		assert entry_chrB.blocks == '1'
		assert entry_chrB.bp_per_block_sum == '50'
		assert entry_chrB.variant_per_block_sum == '2'
		assert entry_chrB.block_n50[:-1] == '0'

		assert entry_all.chromosome == 'ALL'
		assert entry_all.variants == '10'
		assert entry_all.phased == '9'
		assert entry_all.unphased == '1'
		assert entry_all.blocks == '3'
		assert entry_all.bp_per_block_sum == '601'
		assert entry_all.variant_per_block_sum == '9'
		assert entry_all.block_n50[:-1] == '0'

def test_stats2():
	for (filename, ploidy) in [('tests/data/phased3.vcf', 2), ('tests/data/phased3-poly.vcf', 3)]:
		with TemporaryDirectory() as tempdir:
			outtsv = tempdir + '/output.tsv'
			run_stats(ploidy=ploidy, vcf=filename, tsv=outtsv, sample='sample1', chr_lengths='tests/data/chr-lengths.txt')
			lines = [l.split('\t') for l in open(outtsv)]
			assert len(lines) == 4
			Fields = namedtuple('Fields', [f.strip('#\n') for f in lines[0]])
			entry_chrA, entry_chrB, entry_all = [Fields(*l) for l in lines[1:]]

			assert entry_chrA.chromosome == 'chrA'
			assert entry_chrA.variants == '9'
			assert entry_chrA.phased == '4'
			assert entry_chrA.unphased == '5'
			assert entry_chrA.blocks == '1'
			assert entry_chrA.variant_per_block_sum == '4'
			assert entry_chrA.bp_per_block_sum == '350'
			assert entry_chrA.block_n50[:-1] == '0'

			assert entry_chrB.chromosome == 'chrB'
			assert entry_chrB.variants == '4'
			assert entry_chrB.phased == '4'
			assert entry_chrB.unphased == '0'
			assert entry_chrB.blocks == '1'
			assert entry_chrB.variant_per_block_sum == '4'
			assert entry_chrB.bp_per_block_sum == '400'
			assert entry_chrB.block_n50[:-1] == '400'

			assert entry_all.chromosome == 'ALL'
			assert entry_all.variants == '13'
			assert entry_all.phased == '8'
			assert entry_all.unphased == '5'
			assert entry_all.blocks == '2'
			assert entry_all.variant_per_block_sum == '8'
			assert entry_all.bp_per_block_sum == '750'
			assert entry_all.block_n50[:-1] == '350'

def test_stats_bed1():
	for (filename, ploidy) in [('tests/data/phased3.vcf', 2), ('tests/data/phased3-poly.vcf', 3)]:
		with TemporaryDirectory() as tempdir:
			outtsv = tempdir + '/output.tsv'
			run_stats_bed(ploidy=ploidy, vcf=filename, bed='tests/data/phased3.bed', tsv=outtsv, sample='sample2', chr_lengths='tests/data/chr-lengths.txt')
			lines = [l.split('\t') for l in open(outtsv)]
			assert len(lines) == 7
			Fields = namedtuple('Fields', [f.strip('#\n') for f in lines[0]])
			entry_interval1, entry_interval2, entry_chrA, entry_interval3, entry_chrB, entry_all = [Fields(*l) for l in lines[1:]]

			assert entry_interval1.chromosome == 'chrA:49-399'
			assert entry_interval1.variants == '4'
			assert entry_interval1.phased == '4'
			assert entry_interval1.unphased == '0'
			assert entry_interval1.blocks == '2'
			assert entry_interval1.variant_per_block_sum == '4'
			assert entry_interval1.bp_per_block_sum == '102'

			assert entry_interval2.chromosome == 'chrA:499-849'
			assert entry_interval2.variants == '4'
			assert entry_interval2.phased == '4'
			assert entry_interval2.unphased == '0'
			assert entry_interval2.blocks == '1'
			assert entry_interval2.variant_per_block_sum == '4'
			assert entry_interval2.bp_per_block_sum == '302'

			assert entry_chrA.chromosome == 'chrA'
			assert entry_chrA.variants == '8'
			assert entry_chrA.phased == '8'
			assert entry_chrA.unphased == '0'
			assert entry_chrA.blocks == '3'
			assert entry_chrA.variant_per_block_sum == '8'
			assert entry_chrA.bp_per_block_sum == '404'

			assert entry_interval3.chromosome == 'chrB:99-199'
			assert entry_interval3.variants == '2'
			assert entry_interval3.phased == '0'
			assert entry_interval3.unphased == '2'
			assert entry_interval3.blocks == '0'
			assert entry_interval3.variant_per_block_sum == '0'
			assert entry_interval3.bp_per_block_sum == '0'

			assert entry_chrB.chromosome == 'chrB'
			assert entry_chrB.variants == '2'
			assert entry_chrB.phased == '0'
			assert entry_chrB.unphased == '2'
			assert entry_chrB.blocks == '0'
			assert entry_chrB.variant_per_block_sum == '0'
			assert entry_chrB.bp_per_block_sum == '0'
