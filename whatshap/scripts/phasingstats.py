#!/usr/bin/env python3
"""
Read phasing information from a VCF and print out statistics. Optionally,
write out a GTF file that describes the found blocks.

The GTF file can be loaded into IGV. Since a block of phased variants is not
necessarily contiguos, each block is modelled as a "gene" in the GTF file that
can have multiple "exons". Each exon marks a set of adjacent variants which all
are phased. If there are nonadjacent variants that belong to the same
block, there will be multiple "exons" connected with an arrow (in IGV).

The statistics printed at the end exclude blocks of size 1.
"""
"""
TODO
* factor out common code of parse_hp_tags and parse_pipe_notation
"""
import logging
import sys
from collections import Counter, defaultdict
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from time import time
from ..vcf import vcf_sample_reader

__author__ = 'Marcel Martin'

logger = logging.getLogger(__name__)


def frequency_median(frequencies):
	"""
	Given a dictionary of frequencies, return the median.
	If the total no. of values is odd, the left of both
	middle values is returned.
	"""
	m = 0  # partial sum
	middle = (1 + sum(frequencies.values())) // 2
	for length in sorted(frequencies):
		m += frequencies[length]
		if m >= middle:
			return length
	# never reached
	assert False


class GtfWriter:
	def __init__(self, file):
		self._file = file

	def write(self, chromosome, start, stop, name):
		"""
		Write a feature to the GTF. start is 0-based.
		"""
		assert start < stop
		print(chromosome, 'Phasing', 'exon', start + 1, stop, '.', '+', '.',
			'gene_id "{}"; transcript_id "{}.1";'.format(name, name),
			sep='\t', file=self._file)


def parse_pipe_notation(path):
	assert False, "This has not been updated to print a GTF"
	n = 0
	block_start = None
	singletons = 0
	block_lengths = []  # for statistics
	for record in vcf.Reader(filename=path):
		if not record.is_snp:
			continue
		if len(record.ALT) > 1:
			logger.warning("Reading VCFs with multiple ALTs not implemented")
			continue
		call = record.samples[0]
		assert len(record.samples) == 1
		if not call.is_het:
			continue

		if False:  # for debugging
			print(
				'{:10d}'.format(record.start),
				'alleles:', record.alleles,
				call.gt_alleles,
				'phased:', int(call.phased),
				'het:', int(call.is_het),
				'bases:', call.gt_bases
			)
		# Look for blocks of phased variants. According to the VCF specification,
		# the first variant in such blocks is marked as 'unphased'. Since it
		# needs to be considered part of the phased block, we keep track of the
		# position of the previous variant in prev_start.
		if call.phased:
			assert n > 0, "First variant in VCF is phased. This shouldn't be."
			# block_start is None when we're outside a block
			if block_start is None:
				# start a new block, but start it at the position of the
				# previous variant
				block_start = prev_start
				n_variants = 2  # this one and the one at prev_start
			else:
				n_variants += 1
		else:
			if block_start is not None:
				# end of block
				assert prev_start > block_start
				block_length = prev_start - block_start + 1
				print(record.CHROM, block_start, prev_start + 1, 'Len={},N={}'.format(block_length, n_variants), sep='\t')
				block_lengths.append(block_length)
				block_start = None
			else:
				singletons += 1

		prev_start = record.start
		n += 1
	if block_start is not None:
		assert prev_start > block_start
		block_length = prev_start - block_start + 1
		print(record.CHROM, block_start, prev_start + 1, 'Len={},N={}'.format(block_length, n_variants), sep='\t')
		block_lengths.append(block_length)

	# print statistics
	def printe(*args, **kwargs):
		kwargs['file'] = sys.stderr
		print(*args, **kwargs)

	printe('considered variants:', n)
	printe('blocks:', len(block_lengths))
	if blocks:
		printe('longest block:', max(block_lengths))
		printe('shortest block:', min(block_lengths))
		printe('average block: {:.2f}'.format(sum(block_lengths) / len(block_lengths)))
		printe('singletons (unphased variants not within a block):', singletons)


def parse_hp_tags(path, sample, gtfwriter=None):
	"""
	sample -- None means: use the first sample
	"""
	n_phased = 0
	n_records = 0
	blocks = Counter()  # maps (chromosome, block_name) tuples to variant counts
	block_lengths = defaultdict(int)
	prev_phased = False
	prev_block_name = None
	prev_record = None
	last_time = time()
	for sample, record, call in vcf_sample_reader(path, sample):
		n_records += 1
		assert not (not call.is_het and hasattr(call.data, 'HP')), "HP tag for homozygous variant found"

		if not hasattr(call.data, 'HP'):
			block_name = None
		else:
			assert len(call.data.HP) == 2
			# call.data.HP is something like: ['550267-1', '550267-2']
			block_name = call.data.HP[0].split('-')[0]
			blocks[(record.CHROM, block_name)] += 1
			n_phased += 1

		if prev_block_name != block_name or (prev_record is not None and prev_record.CHROM != record.CHROM):
			# some type of transition is occurring here
			if prev_block_name is not None:
				# phased block has ended at previous variant
				if gtfwriter:
					gtfwriter.write(prev_record.CHROM, block_start, prev_record.start + 1, prev_block_name)
				block_lengths[(prev_record.CHROM, prev_block_name)] += prev_record.start + 1 - block_start
			if block_name is not None:
				# phased block starts
				block_start = record.start

		if n_records % 10000 == 0 and time() - last_time > 20:
			logger.info('%s records processed', n_records)
			last_time = time()
		prev_block_name = block_name
		prev_record = record
	if prev_block_name is not None:
		if gtf:
			gtf.write(prev_record.CHROM, block_start, prev_record.start + 1, prev_block_name)
		block_lengths[(prev_record.CHROM, prev_block_name)] += prev_record.start + 1 - block_start
	logger.info('%s records processed', n_records)

	# Print statistics
	block_sizes = sorted(blocks.values())
	n_singletons = sum(1 for size in block_sizes if size == 1)
	block_sizes = [ size for size in block_sizes if size > 1 ]
	block_lengths = sorted(l for l in block_lengths.values() if l > 1)

	WIDTH = 21
	print('Phasing statistics for', path)
	print('Variants in VCF:'.rjust(WIDTH), '{:8d}'.format(n_records))
	print('Phased:'.rjust(WIDTH), '{:8d}'.format(sum(block_sizes)))
	print('Unphased:'.rjust(WIDTH), '{:8d}'.format(n_records - n_phased), '(not considered below)')
	print('Singletons:'.rjust(WIDTH), '{:8d}'.format(n_singletons), '(not considered below)')
	print('Blocks:'.rjust(WIDTH), '{:8d}'.format(len(block_sizes)))
	if block_sizes:
		print()
		print('Block sizes (no. of variants)')
		print('Median block size:'.rjust(WIDTH), '{:8d}    variants'.format(
			block_sizes[len(block_sizes) // 2]))
		print('Average block size:'.rjust(WIDTH), '{:11.2f} variants'.format(
			sum(block_sizes) / len(block_sizes)))
		print('Largest block:'.rjust(WIDTH), '{:8d}    variants'.format(block_sizes[-1]))
		print('Smallest block:'.rjust(WIDTH), '{:8d}    variants'.format(block_sizes[0]))

	if block_lengths:
		print()
		print('Block lengths (basepairs)')
		print('Sum of lengths:'.rjust(WIDTH), '{:8d}    bp'.format(sum(block_lengths)))
		print('Median block length:'.rjust(WIDTH), '{:8d}    bp'.format(
			block_lengths[len(block_lengths) // 2]))
		print('Average block length:'.rjust(WIDTH), '{:11.2f} bp'.format(
			sum(block_lengths) / len(block_lengths)))
		print('Longest block:'.rjust(WIDTH), '{:8d}    bp'.format(block_lengths[-1]))
		print('Shortest block:'.rjust(WIDTH), '{:8d}    bp'.format(block_lengths[0]))


def main():
	logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')

	parser = ArgumentParser(prog='phasingstats', description=__doc__,
		formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('--gtf', default=None,
		help='Write phased blocks to GTF file.')
	parser.add_argument('--sample', metavar='SAMPLE', default=None,
		help='Name of the sample to process. If not given, use first sample '
			'found in VCF.')
	parser.add_argument('vcf', help='VCF file')
	args = parser.parse_args()

	gtfwriter = None
	if args.gtf:
		gtf_file = open(args.gtf, 'wt')
		gtfwriter = GtfWriter(gtf_file)
	#if args.pipeslash:
		#parse_pipe_notation(args.vcf)
	#else:
	parse_hp_tags(args.vcf, args.sample, gtfwriter)
	if gtfwriter:
		gtf_file.close()


if __name__ == '__main__':
	main()
