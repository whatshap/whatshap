"""
Phasing.

* ReadVariant and ReadVariantList don’t really belong here.
TODO: Once core reads are used everywhere, this file can be removed
"""
from collections import namedtuple
import logging
from .core import Read, ReadSet, DPTable

logger = logging.getLogger(__name__)

# List of variants that belong to a single read.
# The variants attribute is a list of ReadVariant objects (see below).
ReadVariantList = namedtuple('ReadVariantList', 'name mapq variants')

# A single variant on a read.
ReadVariant = namedtuple('ReadVariant', 'position base allele quality')


def read_to_coreread(read):
	if type(read.mapq) is int:
		coreread = Read(read.name, read.mapq)
	elif type(read.mapq) is tuple:
		coreread = Read(read.name, min(read.mapq))
	else:
		assert False, 'Strange MAPQ'
	for variant in read.variants:
		# variant is None when there was a "--" in the wif file, which seperates
		# the two parts of a read pair. Not needed in Read/ReadSet.
		if variant is None: continue
		assert variant.allele in ['0', '1'], 'Unknown allele: {}'.format(variant.allele)
		coreread.addVariant(variant.position, variant.base, int(variant.allele), variant.quality)
	return coreread
