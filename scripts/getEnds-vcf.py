#!/usr/bin/env python3
"""
gets the heterozygous snp positions from a vcf file, and then
gathers those snps that coincide with each read end (one read end's
set of positions per line); and also splits ends into their
respective read groups -- Murray Patterson

tailored toward processing GoNL vcf files -- Alex Schoenhuth

output is sorted by read name
"""
from __future__ import print_function
import sys
from collections import namedtuple
try:
	from sqt import HelpfulArgumentParser as ArgumentParser
except:
	from argparse import ArgumentParser
import pysam
import vcf
import gzip


def parse_vcf(path, chromosome, sample_names):
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
			print("reading VCFs with multiple ALTs not correctly implemented", file=sys.stderr)
			continue

		if indices is None:
			indices = [ (i, call.sample) for i, call in enumerate(record.samples) if call.sample in sample_names ]
			if len(indices) == 0:
				print("Error: none of the sample names found in vcf", file=sys.stderr)
				sys.exit(1)
			else:
				outstring = "Found samples "
				for indtup in indices:
					outstring += "%s " % (indtup[1])
				outstring += "in columns "
				for indtup in indices:
					outstring += "%d " % (indtup[0])
				print(outstring, file=sys.stderr)

		het = False
		for index in [x[0] for x in indices]:
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
			print("not a heterozygous SNP for any of the samples, SNP %s" % (record.POS), file=sys.stderr)
			#% (individual, tk[1]), file=sys.stderr)
			continue
		else: # found a heterozygous snp for the individual

			# TODO deal with len(ALT) > 1
			snp_info = [record.POS, record.REF, record.ALT[0]]
			for index in [x[0] for x in indices]:
				"""
				TODO
				v = tk[index].split(':')[0]
				if v in ('.|1', '0/1'): v = '0|1' # just to disambiguate what
				elif v in ('1|.', '1/0'): v = '1|0' # was mentioned above
				snp_info.append(v)
				"""
				snp_info.append('XXX')
			variants.append(snp_info)

	return variants


# list of variants that belong to a single read
ReadVariantList = namedtuple('ReadVariantList', 'name mapq variants')

# a single variant on a read
ReadVariant = namedtuple('ReadVariant', 'position base allele quality')


def read_bam(path, chromosome, variants, mapq_threshold=20):
	# NOTE: we assume that there are only M,I,D,S (no N,H,P,=,X) in any
	# CIGAR alignment of the bam file

	# first we get some header info, etc.
	samfile = pysam.Samfile(path, "rb")

	target_tid = samfile.gettid(chromosome)
	if target_tid < 0:
		print('ERROR: chromosome unknown in BAM file', file=sys.stderr)
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

		# only reads with a nonempty cigar string (i.e., mapped) are considered
		cigar = read.cigar
		if not cigar:
			continue
		#f = rgF[rgMap[read.opt('RG')]]
		# convert from BAM zero-based coords to 1-based
		pos = int(read.pos) + 1

		# since reads are ordered by position, we need not consider
		# positions that are too small
		while i < len(variants) and variants[i][0] < pos:
			i += 1

		c = 0  # hit count
		j = i  # another index into variants
		p = pos
		s = 0  # absolute index into the read string [0..len(read)]
		# assuming that CIGAR contains only M,I,D,S
		read_variants = []
		#print('Processing read', fl, file=sys.stderr)
		for cigar_op, length in cigar:
			#print('  cigar:', cigar_op, length, file=sys.stderr)
			if cigar_op == 0:  # MATCH/MISMATCH # we're in a matching subregion
				s_next = s + length
				p_next = p + length
				r = p + length  # size of this subregion
				# skip over all SNPs that come before this region
				while j < len(variants) and variants[j][0] < p:
					j += 1
				# iterate over all positions in this subregion and
				# check whether any of them coincide with one of the SNPs ('hit')
				while j < len(variants) and p < r:
					if variants[j][0] == p: # we have a hit
						base = read.seq[s:s+1].decode()
						if base == variants[j][1]:
							al = '0'  # REF allele
						elif base == variants[j][2]:
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
			elif cigar_op == 1 : # an insertion
				s += length
			elif cigar_op == 2 : # a deletion
				p += length
			elif cigar_op == 4 : # soft clipping
				s += length
			elif cigar_op == 5 : # hard clipping
				pass
			else:
				print("error: invalid cigar operation:", cigar_op)
				sys.exit(1)

		#fl += " # " + str(c) + " " + str(read.mapq) + " " + "NA"
		if c > 0:
			rvl = ReadVariantList(name=read.qname, mapq=read.mapq, variants=read_variants)
			result.append(rvl)
	return result


"""
def parse_line(line):
	t = line.split()
	name = t[0]
	count = 0
	mapq = 0
	is_unique = 0
	snps = {}
	for i in range(len(t)):
		if t[i] == ":":
			snps[t[i+1]] = [t[i+2], t[i+3], t[i+4]]
		if t[i] == "#":
			count = int(t[i+1])
			mapq = int(t[i+2])  # mapping quality
			is_unique = t[i+3]  # unique flag ('U' is for unique, 'R' is for repetitive, adopted from BWA XT tag)
	assert count == len(snps)
	return name, count, mapq, is_unique, snps
"""

def merge_ends_and_print_result(variants):
	if len(sys.argv) < 2 :
		print("usage : " + str(sys.argv[0]) + " endsFile")
		sys.exit(0)

	f = open(sys.argv[1],"r")

	# get first end
	e = f.readline()
	if not e :
		print("file is empty")
		sys.exit(0)

	# snps maps a position to a list [base, allele, quality]
	# parse the first line
	name, count, mapq, is_unique, snps = parse_line(e)

	# parse the remaining lines
	while True:
		ep = f.readline() # get second end
		if not ep: # no ep: end e is unpaired
			# seems we are at EOF
			if count > 1: # so simply print end e
				for p in sorted(snps.keys()) : # careful: the default is lists in no particular order, but we want snps to be ordered on their fragment
					print(p + " " + snps[p][0] + " " + snps[p][1] + " " + snps[p][2] + " : ", end='')
				print("# " + str(mapq) + " : " + is_unique)
			break

		# everything with the 'p' suffix is from the second (paired) read
		np, cp, mp, up, sp = parse_line(ep)

		if name == np: # end e pairs up with end ep
			# output merged pair (a read)
			uup = 0
			if count + cp > 1:
				for p in sorted(snps.keys()):
					print(p + " " + snps[p][0] + " " + snps[p][1] + " " + snps[p][2] + " : ", end='')
				print("-- : ", end='') # add a symbol for gap in paired-end reads
				for p in sorted(sp.keys()) :
					print(p + " " + sp[p][0] + " " + sp[p][1] + " " + sp[p][2] + " : ", end='')
				uup = "%s %s" % (is_unique,up)
	#            if is_unique == up == 'U': # uniquely mapped if both ends are
	#                uup = 'U'
	#            else:
	#                uup = 'R'
				print("# " + str(mapq) + " " + str(mp) + " : " + uup) # old: str((mapq+mp)/2.0) + " " + uup
				# note: replace avg of mapq's and display both

			# get new end for next iter
			e = f.readline()
			if not e:
				break

			name, count, mapq, is_unique, snps = parse_line(e)
		else:
			if count > 1: # simply print end
				for p in sorted(snps.keys()) :
					print(p + " " + snps[p][0] + " " + snps[p][1] + " " + snps[p][2] + " : ", end='')
				print("# " + str(mapq) + " : " + is_unique)
			else:
				pass
				#print('not printing', snps, file=sys.stderr)
			e = ep # and use ep for end of next iter
			name = np
			count = cp
			mapq = mp
			snps = sp
			is_unique = up



def print_wif(reads):
	for read in reads:
		print(read.name, end='')
		for variant in read.variants:
			print(' : {position} {base} {allele} {quality}'.format(**vars(variant)), end='')
		print(" # {} {} NA".format(len(read.variants), read.mapq))

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
	parser = ArgumentParser(description=__doc__)
	parser.add_argument('bam', metavar='BAM', help='BAM file')
	parser.add_argument('vcf', metavar='VCF', help='VCF file')
	parser.add_argument('chromosome', help='chromosome to consider')
	parser.add_argument('samples', metavar='sample', nargs='+', help='name(s) of the samples to consider')
	args = parser.parse_args()

	variants = parse_vcf(args.vcf, args.chromosome, args.samples)

	print('Read %d SNPs on chromosome %s' % (len(variants), args.chromosome), file=sys.stderr)

	reads_with_variants = read_bam(args.bam, args.chromosome, variants)

	# sort by read name
	reads_with_variants.sort(key=lambda r: r.name)
	#print_wif(reads_with_variants)
	merge_ends_and_print_result(reads_with_variants)


if __name__ == '__main__':
	main()
