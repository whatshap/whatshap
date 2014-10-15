#!/usr/bin/env python3

import sys
#from bitarray import bitarray

def read_wif(filename):
	'''Returns an iterator that returns lists ([(pos,nucleotide,0/1,quality),..], suffix, original_line)'''
	skipped_reads = 0
	total_reads = 0
	for line in (s.strip() for s in open(filename)):
		total_reads += 1
		fields = [x.strip() for x in line.split(':')]
		assert len(fields) > 2
		assert fields[-2].startswith('#')
		suffix = fields[-2:]
		fields = fields[:-2]
		read = []
		skip_read = False
		last_pos = -1
		for field in fields:
			if field == '--': continue
			tokens = field.split()
			assert len(tokens) == 4
			if tokens[2] == 'E':
				skip_read = True
				break
			pos, nucleotide, bit, quality = int(tokens[0]), tokens[1], tokens[2], int(tokens[3])
			assert nucleotide in ['A', 'C', 'G', 'T', '0', '1', '-', 'X']
			if not last_pos < pos:
				skip_read = True
				break
			read.append((pos-1, nucleotide, bit, quality))
			last_pos = pos
		if skip_read:
			skipped_reads += 1
			continue
		yield read, suffix, line
	print('read_wif(%s): skipped'%filename, skipped_reads, 'out of', total_reads, 'reads', file=sys.stderr)

def wif_to_position_list(filename):
	'''Read wif file and return sorted list of positions present in that file.'''
	position_set = set()
	for read, suffix, line in read_wif(filename):
		for pos, nucleotide, bit, quality in read:
			position_set.add(pos)
	position_list = list(position_set)
	position_list.sort()
	return position_list
