#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division
from optparse import OptionParser, OptionGroup
import sys
import os
from wifreader import read_wif, wif_to_position_list

__author__ = "Tobias Marschall"

usage = """%prog [options] <input.wif> <output-prefix>

Reads <input.wif>, assigns reads to slices (based on the order of input),
and outputs files <output-prefix>.<slice-nr>.wif, for as many slices as necessary."""


class CoverageMonitor:
	'''TODO: This is a most simple, naive implementation. Could do this smarter.'''
	def __init__(self, length):
		self.coverage = [0] * length
	def max_coverage_in_range(self, begin, end):
		return max(self.coverage[begin:end])
	def add_read(self, begin, end):
		for i in xrange(begin,end):
			self.coverage[i] += 1

def position_set(read_list):
	result = set()
	for read, suffix, line in read_list:
		for pos, nucleotide, bit, quality in read:
			result.add(pos)
	return result

def read_cmp(read_record1, read_record2):
	read1, suffix1, line1 = read_record1
	read2, suffix2, line2 = read_record2
	return cmp(read1[0][0], read2[0][0])

def main():
	parser = OptionParser(usage=usage)
	parser.add_option("-H", action="store", dest="slice_height", default=15, type=int,
			help='Maximal height (i.e. coverage) of each slice (default: 15).')

	(options, args) = parser.parse_args()
	if len(args) != 2:
		parser.print_help()
		sys.exit(1)

	input_filename = args[0]
	output_prefix = args[1]
	
	position_list = wif_to_position_list(input_filename)
	print('Found %d SNP positions'%len(position_list), file=sys.stderr)
	
	# dictionary to map SNP position to its index
	position_to_index = dict((position,index) for index,position in enumerate(position_list))
	
	# List of slices, start with one empty slice ...
	slices = [[]]
	# ... and the corresponding coverages along each slice
	slice_coverages = [CoverageMonitor(len(position_list))]
	skipped_reads = 0
	accessible_positions = set()
	for read, suffix, line in read_wif(input_filename):
		# Skip reads that cover only one SNP
		if len(read) < 2:
			skipped_reads += 1
			continue
		for pos, nucleotide, bit, quality in read:
			accessible_positions.add(pos)
		begin = position_to_index[read[0][0]]
		end = position_to_index[read[-1][0]] + 1
		slice_id = 0
		while True:
			# Does current read fit into this slice?
			if slice_coverages[slice_id].max_coverage_in_range(begin, end) < options.slice_height:
				slice_coverages[slice_id].add_read(begin,end)
				slices[slice_id].append((read, suffix, line))
				break
			else:
				slice_id += 1
				# do we have to create a new slice?
				if slice_id == len(slices):
					slices.append([])
					slice_coverages.append(CoverageMonitor(len(position_list)))
	print('Skipped %d reads that only covered one SNP ...'%skipped_reads, file=sys.stderr)
	unphasable_snps = len(position_list) - len(accessible_positions)
	print('... %d out of %d SNP positions (%f%%) where only covered by such reads and are thus unphasable'%(unphasable_snps, len(position_list), unphasable_snps*100.0/len(position_list)), file=sys.stderr)
	# sort slices
	for read_list in slices:
		read_list.sort(cmp=read_cmp)
	for slice_id, read_list in enumerate(slices):
		positions_covered = len(position_set(read_list))
		print('Slice %d contains %d reads and covers %d of %d SNP positions (%f%%)'%(slice_id, len(read_list), positions_covered, len(position_list), positions_covered*100.0/len(position_list)), file=sys.stderr)
		slice_file = open('{0}.{1:02}.wif'.format(output_prefix, slice_id), 'w')
		for read, suffix, line in read_list:
			print(line, file=slice_file)
		slice_file.close()

if __name__ == '__main__':
	sys.exit(main())
