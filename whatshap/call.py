"""
Generate candidate SNP positions from samtools mpileup output.
"""

import pysam
import sys
import re
import pyfaidx
from collections import defaultdict
import datetime
#from whatshap.args import HelpfulArgumentParser as ArgumentParser
import logging

logger = logging.getLogger(__name__)

def add_arguments(parser):
	add = parser.add_argument
	add('--minabs', metavar='MIN_ABS', default=3, type=int,
		help='Minimum absolute ALT depth to call a SNP (default: %(default)s).')
	add('--minrel', metavar='MIN_REL', default=0.4, type=float,
		help='Minimum relative ALT depth to call a SNP (default: %(default)s).')
	add('--multi-allelics', default=False, action='store_true',
		help='Also output multi-allelic sites, if not given only the best ALT allele is reported (if unique).')
	add('--sample', metavar='SAMPLE', default=None, 
		help='Put this sample column into VCF (default: output sites-only VCF).')
	add('ref', metavar='REF', help='FASTA with reference genome')
	add('bam', metavar='BAM', help='BAM file')


def validate(args, parser):
	pass


def main(args):
	fasta = pyfaidx.Fasta(args.ref, as_raw=True)

	print('##fileformat=VCFv4.2')
	print('##fileDate={}'.format(datetime.datetime.now().strftime('%Y%m%d')))
	print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
	print('##FILTER=<ID=PASS,Description="All filters passed">')
	header_columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
	if args.sample is not None:
		header_columns += ['FORMAT', args.sample]
	print(*header_columns, sep='\t')

	re_nucleotide = re.compile('[ACGTNacgtn]')
	re_indel = re.compile('[-\\+]([0-9]+)')
	re_ref = re.compile('[,\\.]')
	re_ignore = re.compile('([\\$\\*]|\\^.)')

	bamfile = pysam.AlignmentFile(args.bam, "rb")
	for pileupcolumn in bamfile.pileup(min_mapping_quality=20, min_base_quality=5):
		try:
			pileup = ''.join(pileupcolumn.get_query_sequences(mark_matches=True, mark_ends=True, add_indels=True))
		except AssertionError:
			# TODO: seems to happen when there are large insertions/deletions
			continue
		position = pileupcolumn.reference_pos + 1
		chromosome = pileupcolumn.reference_name
		i = 0
		bases = defaultdict(int)
#		print('\t'.join([chromosome, str(position), pileup]))
		#print('digesting pileup', pileup)
		n = 0
		ref = fasta[chromosome][position-1].upper()
		while i < len(pileup):
			m = re_nucleotide.match(pileup[i:])
			if m is not None:
				bases[pileup[i].upper()] += 1
				n += 1
				#print('  found base:', pileup[i])
				i += 1
				continue
			m = re_indel.match(pileup[i:])
			if m is not None:
				l = int(m.group(1))
				skip = (m.end()-m.start()) + l
				#print('  found indel:', pileup[i:i+skip])
				i += skip
				continue
			m = re_ref.match(pileup[i:])
			if m is not None:
				bases[ref] += 1
				n += 1
				#print('  found REF:', pileup[i:i+skip])
				i += 1
				continue
			m = re_ignore.match(pileup[i:])
			if m is not None:
				skip = (m.end()-m.start())
				#print('  found other things to ignore:', pileup[i:i+skip])
				i += skip
				continue
			assert False
		#print(bases)
		ref_count = bases[ref]
		alts = []
		for base, count in bases.items():
			if base == ref:
				continue
			if (count >= args.minabs) and (count / (count+ref_count) >= args.minrel):
				alts.append( (count, base) )
		alts.sort(reverse=True)
		if len(alts) > 0:
			columns = [chromosome, position, '.',  ref, '.', '.', 'PASS', '.']
			if args.sample is not None:
				columns += ['GT', '.']
			if args.multi_allelics:
				columns[4] = ','.join(base for count, base in alts)
			else:
				# Do we have two equally supported ALT alleles
				if len(alts) > 1 and (alts[0][0] == alts[1][0]):
					columns[4] = 'N'
				else:
					columns[4] = alts[0][1]
			print(*columns, sep='\t')

if __name__ == '__main__':
	main()
