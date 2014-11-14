"""
Phasing.

TODO
* WIF reading and writing is also in this file, but only because it is only
  used here. The plan is to get rid of WIF files altogether.
* ReadVariant and ReadVariantList don’t really belong here.
"""
from collections import namedtuple
import logging
from tempfile import NamedTemporaryFile
import subprocess
from io import StringIO
from .core import Read,ReadSet,DPTable

logger = logging.getLogger(__name__)

# List of variants that belong to a single read.
# The variants attribute is a list of ReadVariant objects (see below).
ReadVariantList = namedtuple('ReadVariantList', 'name mapq variants')

# A single variant on a read.
ReadVariant = namedtuple('ReadVariant', 'position base allele quality')


def read_wif(filename):
	'''Returns an iterator that returns lists ([(pos,nucleotide,0/1,quality),..], suffix, original_line)'''
	skipped_reads = 0
	total_reads = 0
	for line in open(filename):
		line = line.strip()
		total_reads += 1
		fields = [x.strip() for x in line.split(':')]
		assert len(fields) > 2
		assert fields[-2].startswith('#')
		suffix = fields[-2:]
		fields = fields[:-2]
		variants = []
		skip_read = False
		last_pos = -1
		for field in fields:
			if field == '--':
				variants.append(None)
				continue
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
			variants.append(ReadVariant(position=pos-1, base=nucleotide, allele=bit, quality=quality))
			last_pos = pos
		if skip_read:
			skipped_reads += 1
			continue
		yield ReadVariantList(name=None, mapq=None, variants=variants)
	if skipped_reads > 0:
		logger.warn('read_wif(%s): skipped %d out of %d reads.', filename, skipped_reads, total_reads)


def print_wif(reads, file):
	for read in reads:
		paired = False
		for variant in read.variants:
			if variant is None:
				# this marker is used between paired-end reads
				print('-- : ', end='', file=file)
				paired = True
			else:
				print('{position} {base} {allele} {quality} : '.format(
						position=variant.position + 1,
						base=variant.base,
						allele=variant.allele,
						quality=variant.quality),
					end='', file=file)
		if paired:
			print("# {} {} : NA NA".format(read.mapq[0], read.mapq[1]), file=file)
		else:
			print("# {} : NA".format(read.mapq), file=file)

# output columns:
# - read.qname
# - for each SNP that is on this read:
#   - space, colon, space
#   - position
#   - read base at this position
#   - '0' or '1': 0 for reference allele, 1 for alt allele
#   - base quality at this position
# - finally
#   - space, hash, space
#   - no. of SNPs for this read
#   - mapping quality
#   - "NA"

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
		assert variant.allele in ['0','1'], 'Unknown allele: {}'.format(variant.allele)
		coreread.addVariant(variant.position, variant.base, int(variant.allele), variant.quality)
	return coreread

def coreread_to_read(coreread):
	read = ReadVariantList(name=coreread.getName(), mapq=coreread.getMapq(), variants=[])
	for position, base, allele, quality in coreread:
		read.variants.append(ReadVariant(position=position, base=base, allele=str(allele), quality=quality))
	return read

def phase_reads(reads, all_het=False, wif=None, superwif=None):
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
