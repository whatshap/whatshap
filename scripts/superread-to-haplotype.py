#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division
from optparse import OptionParser, OptionGroup
import sys
import os
from wifreader import read_wif, wif_to_position_list, determine_connectivity
#from bitarray import bitarray

__author__ = "Tobias Marschall"

usage = """%prog [options] <super-reads.wif> <positions>

Reads <super-reads.wif> and a list of SNP positions to output and prints
the haplotype (wrt to these positions) as given by the super read. A dash (-)
is printed for positions about which the super read does not make a statement."""

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-O", action="store", dest="original_reads", default=None,
			help='WIF with original reads. If given, a "|" will be put between any two non-gap positions in the haplotype that are not jointly covered by at least on read.')

	(options, args) = parser.parse_args()
	if len(args) != 2:
		parser.print_help()
		sys.exit(1)

	superread_filename = args[0]
	positions_filename = args[1]
	
	position_list = [int(x) for x in open(positions_filename)]
	position_list.sort()
	position_to_index = dict((position,index) for index,position in enumerate(position_list))
	print('Read %d SNP positions from "%s"'%(len(position_list), positions_filename), file=sys.stderr)

	connected = None
	if options.original_reads != None:
		connected = determine_connectivity(options.original_reads, position_list)

	for read, suffix, line in read_wif(superread_filename):
		haplotype = ['-'] * len(position_list)
		for pos, nucleotide, bit, quality in read:
			if pos in position_to_index:
				haplotype[position_to_index[pos]] = str(bit)
			else:
				print('Warning: super read contains unknown SNP position:', pos, file=sys.stderr)
		for i, p, h in zip(range(len(position_list)), position_list, haplotype):
			print(p, h)
			if connected and i < len(position_list) - 1 and not connected[i] and haplotype[i] != '-' and haplotype[i+1] != '-':
				print('---')


		# If information on "SNP deserts" is available, then input "|" symbols to separate unconnected components
		if connected != None:
			for i in xrange(len(connected) - 1, -1, -1):
				if (not connected[i]) and (haplotype[i] != '-') and (haplotype[i+1] != '-'):
					haplotype.insert(i+1, '|')
		print(''.join(haplotype))


if __name__ == '__main__':
	main()
