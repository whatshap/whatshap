#!/usr/bin/env python3
"""
Read a VCF file and output a BED file that shows phased blocks.

Only the first sample is considered.

Only SNPs are considered (indels are treated as if they did not exist).
Homozygous calls are also ignored.

TODO make sure this works when there's more than one chromosome in the VCF
"""
import sys
import argparse

from sqt import HelpfulArgumentParser
import vcf

__author__ = 'Marcel Martin'


def main():
	parser = HelpfulArgumentParser(description=__doc__)
	parser.add_argument('vcf', help='VCF file')
	args = parser.parse_args()

	n = 0
	block_start = None
	singletons = 0
	block_lengths = []  # for statistics
	for record in vcf.Reader(filename=args.vcf):
		if not record.is_snp:
			continue
		assert len(record.ALT) == 1, "reading VCFs with multiple ALTs not implemented"
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
	printe('longest block:', max(block_lengths))
	printe('shortest block:', min(block_lengths))
	printe('average block: {:.2f}'.format(sum(block_lengths) / len(block_lengths)))
	printe('singletons (unphased variants not within a block):', singletons)



if __name__ == '__main__':
	main()
