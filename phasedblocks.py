#!/usr/bin/env python3
import sys
import vcf


for rec in vcf.Reader(filename=sys.argv[1]):
	call = rec.samples[0]
	assert len(rec.samples) == 1
	print(rec.start, call.gt_alleles, 'phased:', int(call.phased), 'bases:', call.gt_bases)
