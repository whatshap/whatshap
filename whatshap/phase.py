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
from .core import Read,ReadSet

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
		# TODO: Why can a variant be None??
		if variant is None: continue
		assert variant.allele in ['0','1'], 'Unknown allele: {}'.format(variant.allele)
		coreread.addVariant(variant.position, variant.base, int(variant.allele), variant.quality)
	return coreread

def phase_reads(reads, all_het=False, wif=None, superwif=None):
	"""
	Phase reads, return superreads. This function runs the phasing algorithm
	by creating a temporary WIF file, running the 'dp' binary and then
	parsing the created "super reads" output file.

	Intermediate files are written to the paths named by wif and superwif. If
	the parameters are None, a name for the temporary files is made up.

	TODO The temporary files are *not* deleted.
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
	read_set.finalize()
	print(read_set)
	if wif is not None:
		wif_path = wif
		wif_file = open(wif_path, 'wt')
	else:
		wif_file = NamedTemporaryFile(mode='wt', suffix='.wif', prefix='whatshap-', delete=False)
		wif_path = wif_file.name
	with wif_file as wif:
		print_wif(reads, wif)
		logger.info('WIF written to %s', wif_path)

	dp_cmdline = ['build/dp'] + (['--all_het'] if all_het else []) + [wif_path]
	logger.info('Running %s', ' '.join(dp_cmdline))
	superread_result = subprocess.check_output(dp_cmdline, shell=False).decode()

	if superwif is not None:
		superwif_path = superwif
		superwif_file = open(superwif_path, 'wt')
	else:
		superwif_file = NamedTemporaryFile(mode='wt', suffix='.superwif', prefix='whatshap-', delete=False)
		superwif_path = superwif_file.name
	with superwif_file as wif:
		wif.write(superread_result)
		logger.info('Super WIF written to %s', superwif_path)

	superreads = read_wif(superwif_path)
	return superreads
