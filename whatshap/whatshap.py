#!/usr/bin/env python3
"""
Read a VCF and a BAM file and write a WIF file to standard output.
The WIF file is ready to be used as input for the 'dp' program.

(old description:
gets the heterozygous snp positions from a vcf file, and then
gathers those snps that coincide with each read end (one read end's
set of positions per line); and also splits ends into their
respective read groups)

Output:
"haplotype string", where

 0: ref allele
 1: alt allele
 -: unphasable: no coverage of read that covers at least 2 SNPs
 X: unphasable: there is coverage, but still not phasable (tie)

TODO
* Perhaps simplify slice_reads() such that it only creates and returns one slice
* allow only one sample name to be passed to parse_vcf
"""
import logging
import sys
import random
import gzip
from collections import namedtuple
from tempfile import NamedTemporaryFile
import subprocess
from io import StringIO
try:
	from sqt import HelpfulArgumentParser as ArgumentParser
except:
	from argparse import ArgumentParser
import pysam
import vcf

__author__ = "Murray Patterson, Alexander Schönhuth, Tobias Marschall, Marcel Martin"

logger = logging.getLogger(__name__)

# List of variants that belong to a single read.
# The variants attribute is a list of ReadVariant objects (see below).
ReadVariantList = namedtuple('ReadVariantList', 'name mapq variants')

# A single variant on a read.
ReadVariant = namedtuple('ReadVariant', 'position base allele quality')

#VcfVariant = namedtuple('VcfVariant', 'position reference_allele alternative_allele')

class VcfVariant:
	def __init__(self, position, reference_allele, alternative_allele):
		self.position = position
		self.reference_allele = reference_allele
		self.alternative_allele = alternative_allele


def parse_vcf(path, chromosome, sample_names):
	"""
	Read a VCF and return a list of variants. Each entry in the returned list is
	a list of VcfVariant objects.

	path -- Path to VCF file
	chromosome -- Chromosome to work on
	sample_names -- A list of sample names (strings). Extract only calls
		belonging to a sample that occurs in this list.
	"""
	variants = []
	index = -1
	indices = None
	for record in vcf.Reader(filename=path):
		if record.CHROM != chromosome:
			# TODO use .fetch to avoid iterating over entire file
			continue
		if not record.is_snp:
			continue
		if len(record.ALT) != 1:
			logger.warn("Reading VCFs with multiple ALTs not implemented")
			continue
		if indices is None:
			indices = [ (i, call.sample) for i, call in enumerate(record.samples) if call.sample in sample_names ]
			if len(indices) == 0:
				logger.error("None of the sample names found in vcf")
				sys.exit(1)
			else:
				outstring = "Found samples "
				for indtup in indices:
					outstring += "%s " % (indtup[1])
				outstring += "in columns "
				for indtup in indices:
					outstring += "%d " % (indtup[0])
				logger.info(outstring)

		het = False
		for index, _ in indices:
			call = record.samples[index]
			if False:
				print(
					'pos {:10d}'.format(record.start),
					'alleles:', record.alleles,
					call.gt_alleles,
					'phased:', int(call.phased),
					'het:', int(call.is_het),
					'bases:', call.gt_bases
				)

			het = het or call.is_het
			# TODO 1/1 and 1/2

			#het = het or tk[index] == '0|1' or tk[index] == '1|0'
			# note: the "." means "0" in the simulated Venter dataset, but
			# it can also mean "don't know" in certain contexts, so you
			# may have to come back to this line of code -- murray
			#het = het or tk[index] == '.|1' or tk[index] == '1|.'
			# note also that we may need also to consider "0/1", "./1"
			# ... etc. in cases where we have unphased data; so keep this
			# in mind also -- murray
		if not het:
			logger.warn("Not a heterozygous SNP for any of the samples, SNP %s", record.start + 1)
			continue
		else:
			# found a heterozygous snp for the individual
			yield VcfVariant(position=record.start, reference_allele=record.REF,
				 alternative_allele=record.ALT[0])
			for index, _ in indices:
				"""
				# TODO what was this originally supposed to do?
				v = tk[index].split(':')[0]
				if v in ('.|1', '0/1'): v = '0|1' # just to disambiguate what
				elif v in ('1|.', '1/0'): v = '1|0' # was mentioned above
				snp_info.append(v)
				"""


def read_bam(path, chromosome, variants, mapq_threshold=20):
	"""
	path -- path to BAM file
	chromosome -- name of chromosome to work on
	variants -- list of Variant objects (obtained from VCF with parse_vcf)

	Return a list of ReadVariantList objects.
	"""
	# NOTE: we assume that there are only M,I,D,S (no N,H,P,=,X) in any
	# CIGAR alignment of the bam file

	# first we get some header info, etc.
	# TODO use a context manager
	samfile = pysam.Samfile(path, "rb")

	target_tid = samfile.gettid(chromosome)
	if target_tid < 0:
		logger.error('Chromosome "%s" unknown in BAM file', chromosome)
		# TODO raise an exception instead?
		sys.exit(1)

	#rgMap = {} # get mapping from each read tech to its group
	#for r in samfile.header['RG'] :
		#rgMap[r['ID']] = r['SM']
	#if(len(rgMap)==0) :
		#print("error : no read groups in BAM header")
		#print("exiting ...")
		#sys.exit(0)

	#rgs = [] # get the (set of) unique read groups
	#for k in rgMap.keys() :
		#rgs.insert(0,rgMap[k])
	#rgs = sorted(set(rgs))

	#rgF = {} # a file for each read group
	#for e in rgs :
		#fName = pf + "-" + str(e) + ".ends"
		#rgF[e] = open(fName,"w");

	# resulting list of ReadVariantList objects
	result = []

	# now we loop through the bam file
	i = 0 # to keep track of position in variants array (which is in order)
	# the assumption is that reads in samfile are ordered by position
	# one can use samfile.fetch() for doing that
	for read in samfile:
		if read.tid != target_tid: continue
		# TODO: handle additional alignments correctly! find out why they are sometimes overlapping/redundant
		if read.flag & 2048 != 0:
			# print('Skipping additional alignment for read ', read.qname)
			continue
		if read.is_secondary:
			continue
		if read.is_unmapped:
			continue
		if read.mapq < mapq_threshold:
			continue
		cigar = read.cigar
		if not cigar:
			continue
		#f = rgF[rgMap[read.opt('RG')]]
		pos = read.pos

		# since reads are ordered by position, we need not consider
		# positions that are too small
		while i < len(variants) and variants[i].position < pos:
			i += 1

		c = 0  # hit count
		j = i  # another index into variants
		p = pos
		s = 0  # absolute index into the read string [0..len(read)]
		read_variants = []
		for cigar_op, length in cigar:
			#print('  cigar:', cigar_op, length, file=sys.stderr)
			if cigar_op == 0:  # MATCH/MISMATCH # we're in a matching subregion
				s_next = s + length
				p_next = p + length
				r = p + length  # size of this subregion
				# skip over all SNPs that come before this region
				while j < len(variants) and variants[j].position < p:
					j += 1
				# iterate over all positions in this subregion and
				# check whether any of them coincide with one of the SNPs ('hit')
				while j < len(variants) and p < r:
					if variants[j].position == p: # we have a hit
						base = read.seq[s:s+1].decode()
						if base == variants[j].reference_allele:
							al = '0'  # REF allele
						elif base == variants[j].alternative_allele:
							al = '1'  # ALT allele
						else:
							al = 'E' # for "error" (keep for stats purposes)
						rv = ReadVariant(position=p, base=base, allele=al, quality=ord(read.qual[s:s+1])-33)
						read_variants.append(rv)
						c += 1
						j += 1
					s += 1 # advance both read and reference
					p += 1
				s = s_next
				p = p_next
			elif cigar_op == 1:  # an insertion
				s += length
			elif cigar_op == 2:  # a deletion
				p += length
			elif cigar_op == 4:  # soft clipping
				s += length
			elif cigar_op == 5:  # hard clipping
				pass
			else:
				logger.error("Invalid CIGAR operation: %d", cigar_op)
				sys.exit(1)
		if c > 0:
			rvl = ReadVariantList(name=read.qname, mapq=read.mapq, variants=read_variants)
			result.append(rvl)
	return result


def grouped_by_name(reads_with_variants):
	"""
	Group an input list of reads into a list of lists where each
	sublist is a slice of the input list that contains reads with the same name.

	Example (only read names are shown here):
	['A', 'A', 'B', 'C', 'C']
	->
	[['A', 'A'], ['B'], ['C', 'C']]
	"""
	result = []
	prev = None
	current = []
	for read in reads_with_variants:
		if prev and read.name != prev.name:
			result.append(current)
			current = []
		current.append(read)
		prev = read
	result.append(current)
	return result


def merge_reads(reads, mincount=2):
	"""
	Merge reads that occur twice (according to their name) into a single read.
	This is relevant for paired-end or mate pair reads.

	The ``variants`` attribute of a merged read contains the variants lists of
	both reads, separated by "None". For a merged read, the ``mapq`` attribute
	is a tuple consisting of the two original mapping qualities.

	Reads that occur only once are returned unchanged (unless mincount
		applies).

	mincount -- If the number of variants on a read occurring once or the
		total number of variants on a paired-end read is lower than this
		value, the read (or read pair) is discarded.

	TODO mincount filtering should be done in a different function.
	"""
	result = []
	for group in grouped_by_name(reads):
		count = sum(len(read.variants) for read in group)
		if count < mincount:
			continue
		if len(group) == 1:
			result.append(group[0])
		elif len(group) == 2:
			merged_variants = group[0].variants + [None] + group[1].variants
			merged_read = ReadVariantList(name=group[0].name, mapq=(group[0].mapq, group[1].mapq), variants=merged_variants)
			result.append(merged_read)
		else:
			assert len(group) <= 2, "More than two reads with the same name found"
	return result


def filter_reads(reads):
	"""
	Return a new list in which reads are omitted that fulfill at least one of
	these conditions:

	- one of the read's variants' alleles is 'E'

	- variant positions are not strictly monotically increasing.
	"""
	result = []
	for read in reads:
		prev_pos = -1
		for variant in read.variants:
			if variant is None:
				continue
			if not prev_pos < variant.position:
				break
			if variant.allele == 'E':
				break
			assert variant.base in 'ACGT01-X', 'variant.base={!r}'.format(variant.base)
			prev_pos = variant.position
		else:
			# executed when no break occurred above
			result.append(read)
	return result


class CoverageMonitor:
	'''TODO: This is a most simple, naive implementation. Could do this smarter.'''
	def __init__(self, length):
		self.coverage = [0] * length

	def max_coverage_in_range(self, begin, end):
		return max(self.coverage[begin:end])

	def add_read(self, begin, end):
		for i in range(begin, end):
			self.coverage[i] += 1


def position_set(reads):
	"""
	Return a set of all variant positions that occur within a list of reads.

	reads -- a list of ReadVariantList objects
	"""
	positions = set()
	for read in reads:
		positions.update(variant.position for variant in read.variants if variant is not None)
	return positions


def slice_reads(reads, max_coverage):
	"""
	TODO document this
	TODO document that fragment (not read) coverage is used

	max_coverage --
	reads -- a list of ReadVariantList objects
	"""
	position_list = sorted(position_set(reads))
	logger.info('Found %d SNP positions', len(position_list))

	# dictionary to map SNP position to its index
	position_to_index = { position: index for index, position in enumerate(position_list) }

	# List of slices, start with one empty slice ...
	slices = [[]]
	# ... and the corresponding coverages along each slice
	slice_coverages = [CoverageMonitor(len(position_list))]
	skipped_reads = 0
	accessible_positions = set()
	for read in reads:
		# Skip reads that cover only one SNP
		if len(read.variants) < 2:
			skipped_reads += 1
			continue
		for variant in read.variants:
			if variant is None:
				continue
			accessible_positions.add(variant.position)
		begin = position_to_index[read.variants[0].position]
		end = position_to_index[read.variants[-1].position] + 1
		slice_id = 0
		while True:
			# Does current read fit into this slice?
			if slice_coverages[slice_id].max_coverage_in_range(begin, end) < max_coverage:
				slice_coverages[slice_id].add_read(begin, end)
				slices[slice_id].append(read)
				break
			else:
				slice_id += 1
				# do we have to create a new slice?
				if slice_id == len(slices):
					slices.append([])
					slice_coverages.append(CoverageMonitor(len(position_list)))
	logger.info('Skipped %d reads that only cover one SNP', skipped_reads)

	unphasable_snps = len(position_list) - len(accessible_positions)
	logger.info('%d out of %d variant positions (%.1d%%) do not have a read '
		'connecting them to another variant and are thus unphasable',
		unphasable_snps, len(position_list),
		100. * unphasable_snps / len(position_list))

	# Sort each slice
	for read_list in slices:
		read_list.sort(key=lambda r: r.variants[0].position)
	# Print stats
	for slice_id, read_list in enumerate(slices):
		positions_covered = len(position_set(read_list))
		logger.info('Slice %d contains %d reads and covers %d of %d SNP positions (%.1f%%)',
			  slice_id, len(read_list), positions_covered, len(position_list),
			  positions_covered * 100.0 / len(position_list))

	return slices


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
		yield read
	if skipped_reads > 0:
		logger.warn('read_wif(%s): skipped %d out of %d reads.', filename, skipped_reads, total_reads)


def determine_connectivity(wif_filename, position_list):
	'''Reads WIF of original reads and return a bitarray where bit i says whether
	positions i and i+1 are jointly covered by a read.'''
	position_to_index = dict((position,index) for index,position in enumerate(position_list))
	#b = bitarray(len(position_list)-1)
	#b.setall(0)
	b = [False]*(len(position_list)-1)
	for read in read_wif(wif_filename):
		try:
			start = position_to_index[read[0][0]]
			end = position_to_index[read[-1][0]]
		except KeyError:
			continue
		for i in range(start,end):
			b[i] = True
	return b


def superread_to_haplotype(superread_path, position_list, original_reads):
	position_list = sorted(position_list)
	position_to_index = dict((position,index) for index, position in enumerate(position_list))
	connected = determine_connectivity(original_reads, position_list)

	for read in read_wif(superread_path):
		haplotype = ['-'] * len(position_list)
		for pos, nucleotide, bit, quality in read:
			if pos in position_to_index:
				haplotype[position_to_index[pos]] = str(bit)
			else:
				logger.warn('Super read contains unknown SNP position: %d', pos)
		for i, (p, h) in enumerate(zip(position_list, haplotype)):
			print(p+1, h)
			if connected and i < len(position_list) - 1 and not connected[i] and haplotype[i] != '-' and haplotype[i+1] != '-':
				print('---')

		# If information on "SNP deserts" is available, then input "|" symbols to separate unconnected components
		for i in range(len(connected) - 1, -1, -1):
			if (not connected[i]) and (haplotype[i] != '-') and (haplotype[i+1] != '-'):
				haplotype.insert(i+1, '|')
		print(''.join(haplotype))


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


def main():
	logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
	parser = ArgumentParser(description=__doc__)
	parser.add_argument('--max-coverage', '-H', metavar='MAXCOV', default=15, type=int,
		help='Reduce coverage to at most MAXCOV (default: %(default)s).')
	parser.add_argument('--seed', default=123, type=int, help='Random seed (default: %(default)s)')
	parser.add_argument('--all-het', action='store_true', default=False,
		help='Assume all positions to be heterozygous (that is, fully trust SNP calls).')
	parser.add_argument('--wif', metavar='WIF', default=None, help='Write intermediate WIF file')
	parser.add_argument('--superwif', metavar='SUPERWIF', default=None,
		help='Write intermediate SUPERWIF file')
	parser.add_argument('--resume-wif', metavar='WIF', default=None,
		help='Do not compute WIF, but read it from WIF.')
	parser.add_argument('--resume-superwif', metavar='SUPERWIF', default=None,
		help='Do not compute super WIF, but read it from SUPERWIF.')
	parser.add_argument('bam', metavar='BAM', help='BAM file')
	parser.add_argument('vcf', metavar='VCF', help='VCF file')
	parser.add_argument('chromosome', help='Chromosome to work on')
	parser.add_argument('samples', metavar='SAMPLE', nargs='+', help='Name(s) of the samples to consider')
	args = parser.parse_args()

	if bool(args.resume_superwif) != bool(args.resume_wif):
		parser.error('When resuming, both --resume-wif and --resume-superwif '
			'are required.')
	variants = list(parse_vcf(args.vcf, args.chromosome, args.samples))
	logger.info('Read %d SNPs on chromosome %s', len(variants), args.chromosome)

	if args.resume_wif is None:
		reads_with_variants = read_bam(args.bam, args.chromosome, variants)
		reads_with_variants.sort(key=lambda read: read.name)
		reads = merge_reads(reads_with_variants)

		# sort by position of first variant
		#reads.sort(key=lambda read: read.variants[0].position)

		random.seed(args.seed)
		random.shuffle(reads)
		unfiltered_length = len(reads)
		reads = filter_reads(reads)
		logger.info('Filtered reads: %d', unfiltered_length - len(reads))
		reads = slice_reads(reads, args.max_coverage)[0]

		if args.wif is not None:
			wif_path = args.wif
			wif_file = open(wif_path, 'wt')
		else:
			wif_file = NamedTemporaryFile(mode='wt', suffix='.wif', prefix='whatshap-', delete=False)
			wif_path = wif_file.name
		with wif_file as wif:
			print_wif(reads, wif)
			logger.info('WIF written to %s', wif_path)

		dp_cmdline = ['build/dp'] + (['--all_het'] if args.all_het else []) + [wif_path]
		logger.info('Running %s', ' '.join(dp_cmdline))
		superread_result = subprocess.check_output(dp_cmdline, shell=False).decode()

		if args.superwif is not None:
			superwif_path = args.superwif
			superwif_file = open(superwif_path, 'wt')
		else:
			superwif_file = NamedTemporaryFile(mode='wt', suffix='.superwif', prefix='whatshap-', delete=False)
			superwif_path = superwif_file.name
		with superwif_file as wif:
			wif.write(superread_result)
			logger.info('Super WIF written to %s', superwif_path)
	else:
		superwif_path = args.resume_superwif
		wif_path = args.resume_wif

	positions = [ variant.position for variant in variants ]
	superread_to_haplotype(superwif_path, positions, wif_path)


if __name__ == '__main__':
	main()
