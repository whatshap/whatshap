#!/usr/bin/env python3

from optparse import OptionParser, OptionGroup
import sys
import os
from wifreader import read_wif, wif_to_position_list, determine_connectivity
#from bitarray import bitarray
import vcf


__author__ = "Tobias Marschall"

usage = """%prog [options] <super-reads.wif> <vcf> <chromosome>

Reads <super-reads.wif> and a list of SNP positions to output and prints
the haplotype (wrt to these positions) as given by the super read. A dash (-)
is printed for positions about which the super read does not make a statement."""


allowed_dna_chars = frozenset('ACGTNacgtn')


def read_positions(vcfpath, chromosome):
	positions = []
	for record in vcf.Reader(filename=vcfpath):
		if record.CHROM != chromosome:
			continue
		if not record.is_snp:
			continue
		assert len(record.samples) == 1
		call = record.samples[0]
		if not call.is_het:
			continue
		if len(record.ALT) > 1:
			continue
		if not (set(record.REF) <= allowed_dna_chars):
			continue
		if not (set(record.ALT[0].sequence) <= allowed_dna_chars):
			continue
		positions.append(record.POS)
	return positions


def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-O", action="store", dest="original_reads", default=None,
			help='WIF with original reads. If given, a "|" will be put between any two non-gap positions in the haplotype that are not jointly covered by at least on read.')

	(options, args) = parser.parse_args()
	if len(args) != 1:
		parser.print_help()
		sys.exit(1)

	superread_filename = args[0]
	vcf_filename = args[1]
	chromosome = args[2]

	position_list = read_positions(vcf_filename, chromosome)
	position_list.sort()
	#position_list = list(map(int, '10 20 30 40 50 60 90 100 110 140 150 160 170 180 190'.split()))

	position_to_index = dict((position,index) for index,position in enumerate(position_list))
	print('Read %d SNP positions from "%s"' % (len(position_list), vcf_filename), file=sys.stderr)

	connected = None
	if options.original_reads is not None:
		connected = determine_connectivity(options.original_reads, position_list)

	for read, suffix, line in read_wif(superread_filename):
		haplotype = ['-'] * len(position_list)
		for pos, nucleotide, bit, quality in read:
			if pos in position_to_index:
				haplotype[position_to_index[pos]] = str(bit)
			else:
				print('Warning: super read contains unknown SNP position:', pos, file=sys.stderr)
		for i, (p, h) in enumerate(zip(position_list, haplotype)):
			print(p, h)
			if connected and i < len(position_list) - 1 and not connected[i] and haplotype[i] != '-' and haplotype[i+1] != '-':
				print('---')


		# If information on "SNP deserts" is available, then input "|" symbols to separate unconnected components
		if connected != None:
			for i in range(len(connected) - 1, -1, -1):
				if (not connected[i]) and (haplotype[i] != '-') and (haplotype[i+1] != '-'):
					haplotype.insert(i+1, '|')
		print(''.join(haplotype))


if __name__ == '__main__':
	main()
