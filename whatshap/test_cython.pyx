import math
import logging
import vcf
import cyvcf2
from collections import defaultdict
from .coverage import CovMonitor
from .graph import ComponentFinder
from .priorityqueue import PriorityQueue

from libcpp.unordered_set cimport unordered_set
from .priorityqueue cimport priority_type, priority_type_ptr, queue_entry_type, PriorityQueue
from core cimport ReadSet
cimport cpp

logger = logging.getLogger(__name__)

cdef compute(int bla, int blub):
	cdef int result = 0
	result = bla + blub
	return result

def testcython(vcf_reader2):
#def testcython(outputname):
#	haps1, haps2 = [], []
#	for record in vcf_reader2:
#		for sample in vcf_reader2.samples:
#			if record.genotype(sample).phased:
#				haps1.append(record.genotype(sample)['GT'][0])			
#	#			haps2.append(record.genotype(sample)['GT'][2])
#	print("hap1: ", len(haps1))
	#vcf = VCF(outputname)
#	variant_list = []	
#	for variant in vcf:
#		variant_list.append(variant.POS)
#	print("variant list length: ", len(variant_list))
	cdef int bla = 2
	cdef int blub = 1
	cdef int res = 0
	res = compute(bla, blub)
	print("res: ", res)