#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division
from optparse import OptionParser, OptionGroup
import sys
import os

__author__ = "Tobias Marschall"

usage = """%prog [options] <chromosome> <input.vcf>

Prints all heterozygous SNP positions on a given chromosome."""

allowed_dna_chars = set(['A','C','G','T','N','a','c','g','t','n'])

def valid_dna_string(s):
	chars = set(c for c in s)
	return chars.issubset(allowed_dna_chars)

def main():
	parser = OptionParser(usage=usage)

	(options, args) = parser.parse_args()
	if (len(args) != 2):
		parser.print_help()
		sys.exit(1)

	target_chromosome = args[0]
	vcf_filename = args[1]

	for line in (s.strip() for s in open(vcf_filename)):
		if line.startswith('#'): continue
		fields = line.split()
		chromosome = fields[0]
		position = int(fields[1])
		ref = fields[3]
		alt = fields[4]
		genotype = fields[9].split(':')[0]
		if chromosome != target_chromosome: continue
		if len(ref) != 1: continue
		if len(alt) != 1: continue
		if not valid_dna_string(ref): continue
		if not valid_dna_string(alt): continue
		if not genotype in ['1/0', '0/1', '1|0', '0|1', '1|.', '.|1', '1/.', './1']: continue
		print(position)

if __name__ == '__main__':
	sys.exit(main())
