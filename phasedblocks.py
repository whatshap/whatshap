#!/usr/bin/env python3
"""
Read a VCF file and output a BED file that shows phased blocks.

Only the first sample is considered.

Only SNPs are considered (indels are treated as if they were not existing).
Homozygous calls are also ignored.
"""
import sys
import argparse

from sqt import HelpfulArgumentParser
import vcf


def get_argument_parser():
	return parser


def main():
	parser = HelpfulArgumentParser(description=__doc__)
	parser.add_argument('vcf', help='VCF file')
	args = parser.parse_args()

	n = 0
	within_black = False
	for record in vcf.Reader(filename=args.vcf):
		if not record.is_snp:
			continue
		assert len(record.ALT) == 1, "reading VCFs with multiple ALTs not implemented"
		call = record.samples[0]
		assert len(record.samples) == 1
		if not call.is_het:
			continue

		print('{:10d}'.format(record.start), 'alleles:', record.alleles, call.gt_alleles, 'phased:', int(call.phased), 'het:', int(call.is_het), 'bases:', call.gt_bases)

		previous_pos = record.start
		n += 1
		if n == 50:
			break


if __name__ == '__main__':
	main()
