"""
Phasing.

* ReadVariant and ReadVariantList don’t really belong here.
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


def coreread_to_read(coreread):
	read = ReadVariantList(name=coreread.getName(), mapq=coreread.getMapq(), variants=[])
	for position, base, allele, quality in coreread:
		read.variants.append(ReadVariant(position=position, base=base, allele=str(allele), quality=quality))
	return read


def phase_reads(reads, all_het=False):
	"""
	Phase reads, return superreads. This function runs the phasing algorithm
	via the C++ wrapper.
	"""
	if not reads:
		return [
			ReadVariantList(name=None, mapq=None, variants=[]),
			ReadVariantList(name=None, mapq=None, variants=[])
		]
	# Transform given reads into a core.ReadSet
	read_set = ReadSet()
	for read in reads:
		read_set.add(read_to_coreread(read))
	# Finalizing a read set will sort reads, variants within reads and assign unique read IDs.
	read_set.finalize()

	# Run the core algorithm: construct DP table ...
	dp_table = DPTable(read_set, all_het)
	# ... and do the backtrace to get the solution
	superreads = dp_table.getSuperReads()
	
	# Convert corereads back to "regular" reads
	return  [coreread_to_read(superread) for superread in superreads]
