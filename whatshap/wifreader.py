#!/usr/bin/env python3

import sys

def wif_to_position_list(filename):
	'''Read wif file and return sorted list of positions present in that file.'''
	position_set = set()
	for read, suffix, line in read_wif(filename):
		for pos, nucleotide, bit, quality in read:
			position_set.add(pos)
	position_list = list(position_set)
	position_list.sort()
	return position_list
