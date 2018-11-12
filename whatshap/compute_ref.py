# -*- coding: utf-8 -*-
from __future__ import print_function
import sys
import logging
import itertools
from itertools import chain
import vcf
from cyvcf2 import VCF, Writer
from collections import defaultdict
import collections
from collections import deque
from operator import itemgetter
import customcontainer
import gc
import math
from sys import getsizeof, stderr
import numpy as np
try:
    from reprlib import repr
except ImportError:
    pass

def compute_ref(ref_file, target_file, ref_output):
	#create the list of variants to be considered	
	variant_list = []
	vcf = VCF(target_file)	
	for variant in vcf:
		variant_list.append(variant.POS)

	#compute size of ref_matrix and initialize it with zeros
	#width=number of variants between first and last relevant variant
	#height=twice the number of samples, as each sample offers two haplotypes
	
	vcf_ref = VCF(ref_file)
	height = 2*sum(1 for sample in vcf_ref.samples)
	samplenumber = 0
	for sample in vcf_ref.samples:
		samplenumber+=1
	width = 0
	z = 0
	start = 0
	end = 0
	ref_variant_list = []

	for variant in vcf_ref:
		if (variant.POS == variant_list[0]):
			start = z
		if (variant.POS == variant_list[len(variant_list)-1]):
			end = z
		z += 1
		ref_variant_list.append(variant.POS)
	width = sum(1 for variant in ref_variant_list[start:end+1])

	#fill the matrix
	varcounter = 0
	var2 = 0
	strlist = []
	string = ''
	for i in range(0,height):
		strlist.append(string)
	vcf_ref = VCF(ref_file)
	doublelist = []
	for variant in vcf_ref:
		if (variant.POS in variant_list and variant.POS not in doublelist):
			doublelist.append(variant.POS)	
			var2 += 1
			counter = 0
			for base in variant.gt_bases:			
				if (varcounter >= start and varcounter <= end):
					if (base[0] == variant.REF):
						strlist[counter]+="0"
					else:		
						strlist[counter]+="1"
					if (base[2] == variant.REF):
						strlist[counter+1]+="0"
					else:
						strlist[counter+1]+="1"
				counter += 2
		varcounter += 1
	
	output = open(ref_output, 'w')	
	for string in strlist:
		output.write(string+'\n')
	output.close()
	return output

if (__name__ == '__main__'):
	compute_ref(sys.argv[1], sys.argv[2],sys.argv[3])