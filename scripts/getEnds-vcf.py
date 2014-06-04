#!/usr/bin/env python

# gets the heterozygous snp positions from a vcf file, and then
# gathers those snps that coincide with each read end (one read end's
# set of positions per line); and also splits ends into their
# respective read groups -- Murray Patterson

# tailored toward processing GoNL vcf files -- Alex Schoenhuth

from __future__ import print_function
import sys
import pysam
import gzip


def parse_vcf(path, chromosome, individuals):
	# vcf file
	if path.split(".")[len(path.split("."))-1] == "gz" :
		vcfFile = gzip.open(path,"r") # file is gz'd
	else : vcfFile = open(path,"r") # else it is not

	snpPos = []
	index = -1
	indices = []
	for line in vcfFile.readlines() :
		line = line.strip()
		if line[:2] == '##':
			#print("header line", file=sys.stderr)
			continue
		tk = line.split()
		if tk[0] == '#CHROM':
			print("Determining indices", file=sys.stderr)
			for j, item in enumerate(tk):
				if item in individuals:
					indices.append((j,item))
					#break
			if len(indices) == 0:
				print("Error: none of the individuals found in vcf", file=sys.stderr)
				sys.exit(1)
			else:
				outstring = "Found individuals "
				for indtup in indices:
					outstring += "%s " % (indtup[1])
				outstring += "in columns "
				for indtup in indices:
					outstring += "%d " % (indtup[0])
				print(outstring, file=sys.stderr)
			continue
		if tk[0] != chromosome:
			continue
		if len(tk[3]) != 1 or len(tk[4]) != 1: # no snp in vcf
			#print("No Snp line found, %s" % (tk[1]), file=sys.stderr) # just to avoid polluting stderr
			continue
		#print(tk[index], file=sys.stderr)
		het = False
		for index in [x[0] for x in indices]:
			hetinfo = tk[index].split(':')[0]
			het = het or hetinfo in ['0|1', '1|0', '.|1', '1|.', '0/1', '1/0']
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
			print("not a heterozygous snp for any of the individuals, snp %s" % (tk[1]), file=sys.stderr)
			#% (individual, tk[1]), file=sys.stderr)
			continue
		else: # found a heterozygous snp for the individual

			# tk[0]: chrom
			# tk[1]: pos
			# tk[2]: id
			# tk[3]: ref
			# tk[4]: alt
			snp_info = [int(tk[1]), tk[3], tk[4]]
			#print('snpPos:', snpPos, file=sys.stderr)
			for index in [x[0] for x in indices]:
				v = tk[index].split(':')[0]
				if v in ('.|1', '0/1'): v = '0|1' # just to disambiguate what
				elif v in ('1|.', '1/0'): v = '1|0' # was mentioned above
				snp_info.append(v)
			#print(snpPos[i])
			snpPos.append(snp_info)

	return snpPos


def read_bam_and_print_result(path, chromosome, snpPos):
	# bam file

	# first we get some header info, etc.
	samFile = pysam.Samfile(path,"rb")

	target_tid = samFile.gettid(chromosome)
	if target_tid < 0:
		print('Error chromosome unknown in BAM file', file=sys.stderr)
		sys.exit(1)

	#rgMap = {} # get mapping from each read tech to its group
	#for r in samFile.header['RG'] :
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


	lSnps = len(snpPos)

	# now we loop through the bam file
	i = 0 # to keep track of position in snpPos array (which is in order)
	# the assumption is that reads in samFile are ordered by position
	# one can use samFile.fetch() for doing that
	for read in samFile.fetch() :
		if read.tid != target_tid: continue
		# TODO: handle additional alignments correctly! find out why they are sometimes overlapping/redundant
		if read.flag & 2048 != 0:
			#print('Skipping additional alignment for read ', read.qname)
			continue
		if read.is_secondary: continue
		if read.is_unmapped: continue
		if read.mapq < 20: continue

		# only reads with a nonempty cigar string (i.e., mapped) are considered
		cigar = read.cigar
		if not cigar:
			continue
		#f = rgF[rgMap[read.opt('RG')]]
		# convert from BAM zero-based coords to 1-based
		pos = int(read.pos)+1

		# since reads are ordered by position, we need not consider
		# positions that are too small
		while i < lSnps and snpPos[i][0] < pos:
			i += 1

		c = 0 # hit count
		j = i # another index into snpPos
		p = pos
		s = 0 # absolute index into the read string [0..len(read)]
		# assuming that CIGAR contains only M,I,D,S
		fl = read.qname
		#print('Processing read', fl, file=sys.stderr)
		for cigar_op, length in cigar:
			#print('  cigar:', cigar_op, length, file=sys.stderr)
			if cigar_op == 0:  # MATCH/MISMATCH # we're in a matching subregion
				s_next = s + length
				p_next = p + length
				r = p + length  # size of this subregion
				# skip over all SNPs that come before this region
				while j < lSnps and snpPos[j][0] < p:
					j += 1
				# iterate over all positions in this subregion and
				# check whether any of them coincide with one of the SNPs ('hit')
				while j < lSnps and p < r:
					if snpPos[j][0] == p: # we have a hit
						if read.seq[s] == snpPos[j][1]:
							al = '0'  # REF allele
						elif read.seq[s] == snpPos[j][2]:
							al = '1'  # ALT allele
						else:
							al = 'E' # for "error" (keep for stats purposes)
						fl += " : " + str(p) + " " + str(read.seq[s]) + " " + al + " " + str(ord(read.qual[s])-33) # 34 is the phred score offset in bam files
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
		fl += " # " + str(c) + " " + str(read.mapq) + " " + "NA"
		if c > 0: print(fl)

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

	if len(sys.argv) < 5:
		print("usage : " + str(sys.argv[0]) + " [bam file] [.vcf file] [chromosome] [VCF individual(s)]", file=sys.stderr)
		sys.exit(1)

	# NOTE: we assume that there are only M,I,D,S (no N,H,P,=,X) in any
	# CIGAR alignment of the bam file
	bam = sys.argv[1] #bam = "/data1/gonl/chr22.bam"
	vcfName = sys.argv[2] #vcfName = "gonl.release4.chr22.annotated.vcf.gz"
	chromosome = sys.argv[3]
	individuals = sys.argv[4:] #individuals = "A21c" or = "A2a A2b", etc.

	snpPos = parse_vcf(vcfName, chromosome, individuals)

	# snpPos is a list. each entry is a list: [pos, ref, alt, hetinfo]

	lSnps = len(snpPos)
	print('Read %d SNPs on chromosome %s'%(lSnps,chromosome), file=sys.stderr)
	#for i in range(len(snpPos)) :
		#print(snpPos[i])

	#sys.exit(0)

	read_bam_and_print_result(bam, chromosome, snpPos)


if __name__ == '__main__':
	main()
