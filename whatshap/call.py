"""
Call variants
"""
import logging
import pysam
import pyfaidx
from collections import defaultdict, namedtuple, deque
from contextlib import ExitStack
from ._call import enumerate_reference_kmers

logger = logging.getLogger(__name__)


def add_arguments(parser):
	add = parser.add_argument
	#add('--sample', metavar='SAMPLE', default=None, help='Name of the sample '
			#'to process. If not given, use first sample found in VCF.')
	#add('--chromosome', dest='chromosome', metavar='CHROMOSOME', default=None,
		#help='Name of chromosome to process. If not given, all chromosomes in the '
		#'input VCF are considered.')
	add('reference', metavar='FASTA', help='Reference genome FASTA')
	add('bam', metavar='BAM', help='BAM file ')

DNA_TO_HASH = {'A':0, 'C':1, 'G':2, 'T':3}
HASH_TO_DNA = dict((v,k) for k,v in DNA_TO_HASH.items())

BAM_CMATCH = 0     # M
BAM_CINS = 1       # I
BAM_CDEL = 2       # D
BAM_CREF_SKIP = 3  # N
BAM_CSOFT_CLIP = 4 # S
BAM_CHARD_CLIP = 5 # H
BAM_CPAD = 6       # P
BAM_CEQUAL = 7     # =
BAM_CDIFF = 8      # X
BAM_CBACK = 9      # B

def hash_to_dna(h, k):
	l = []
	for i in range(k):
		l.append(HASH_TO_DNA[h & 3])
		h = h >> 2
	return ''.join(l[::-1])


#def enumerate_reference_kmers(reference, k):
	#h = 0
	#mask = (1 << (2*k)) - 1
	#i = 0
	#for i in range(len(reference)):
		## update hash
		#h = ((h << 2) | DNA_TO_HASH[reference[i]]) & mask
		#if i >= k-1:
			#yield (h, i)


def enumerate_kmers(bam_alignment, k):
	'''
	Generates all kmers in the read and yields pairs (kmer_hash, position),
	where kmer_hash is a binary representation of the kmer and postion is
	the reference position this kmer has been aligned to.
	'''
	h = 0
	mask = (1 << (2*k)) - 1
	pos = bam_alignment.pos
	cigar_index = 0
	cigar_op, cigar_length = bam_alignment.cigartuples[cigar_index]
	i = 0
	consecutive = 0
	while i < len(bam_alignment.query):
		# process cigar entries that don't consume a character from the read
		while True:
			if (cigar_op == BAM_CDEL) or (cigar_op == BAM_CREF_SKIP):
				#print(pos, '-->', pos + cigar_length)
				pos += cigar_length
			elif cigar_op == BAM_CSOFT_CLIP:
				#i += cigar_length
				consecutive = 0
			elif (cigar_length == 0) or (cigar_op == BAM_CHARD_CLIP):
				pass
			else:
				break
			cigar_index += 1
			cigar_op, cigar_length = bam_alignment.cigartuples[cigar_index]

		if i >= len(bam_alignment.query):
			break

		# update hash
		#print('i', i, 'pos', pos, bam_alignment.query[i])
		h = ((h << 2) | DNA_TO_HASH[bam_alignment.query[i]]) & mask
		consecutive += 1
		if consecutive >= k:
			yield (h, pos+1)

		# consume one character of read
		assert cigar_length > 0
		if (cigar_op == BAM_CMATCH) or (cigar_op == BAM_CEQUAL) or (cigar_op == BAM_CDIFF):
			cigar_length -= 1
			pos += 1
		elif cigar_op == BAM_CINS:
			cigar_length -= 1
		else:
			assert False, 'Unexpected cigar operation'

		i += 1


class Caller:
	def __init__(self, reference, k):
		self.reference = reference
		self.k = k
		# bam record for each active read
		self.bam_records = deque()
		# k-mer generators for each active read
		self.kmer_generators = deque()
		self.kmer_generators_finished = deque()
		# .. corresponding queue with the latest (kmer, position) for each read
		self.current_kmers = deque()
		self.ref_kmer_generator = enumerate_reference_kmers(str(reference).encode('UTF-8'),k)
		# position that every other operation is relative to
		self.pileup_columns = deque()
		self.ref_kmers = deque()
		kmer, pos = self.ref_kmer_generator.__next__()
		self.pileup_columns.append(defaultdict(int))
		self.ref_kmers.append(kmer)
		self.ref_pos = pos

	def add_read(self, bam_alignment):
		self.bam_records.append(bam_alignment)
		self.kmer_generators.append(enumerate_kmers(bam_alignment,self.k))
		self.kmer_generators_finished.append(False)
		kmer, pos = self.kmer_generators[-1].__next__()
		self.current_kmers.append((kmer,pos))
		pileup_column, ref_kmer = self.get_column(pos)
		pileup_column[kmer] += 1
		#print(self.current_kmers)

	def finish(self):
		pass

	def get_column(self, pos):
		index = pos - self.ref_pos
		while len(self.pileup_columns) <= index:
			kmer, pos = self.ref_kmer_generator.__next__()
			self.ref_kmers.append(kmer)
			self.pileup_columns.append(defaultdict(int))
		return self.pileup_columns[index], self.ref_kmers[index]

	def pop_column(self):
		if len(self.pileup_columns) > 0:
			result = (self.ref_pos, self.ref_kmers.popleft(), self.pileup_columns.popleft())
		else:
			kmer, pos = self.ref_kmer_generator.__next__()
			assert pos == self.ref_pos
			result = (self.ref_pos, kmer, None)
		self.ref_pos += 1
		return result

	def advance_to(self, target_pos):
		'''
		Add all k-mer from all reads up to target_pos to pileup_columns.
		'''
		#print('advance_to', target_pos)
		for i, kmer_generator in enumerate(self.kmer_generators):
			try:
				kmer, pos = self.current_kmers[i]
				while pos < target_pos:
					kmer, pos = kmer_generator.__next__()
					#print('  read', i, 'kmer', hash_to_dna(kmer, self.k), 'pos', pos)
					pileup_column, ref_kmer = self.get_column(pos)
					pileup_column[kmer] += 1
				self.current_kmers[i] = (kmer,pos)
			except StopIteration:
				self.kmer_generators_finished[i] = True
				pass
		while (len(self.kmer_generators) > 0) and self.kmer_generators_finished[0]:
			#print('  popping read')
			self.current_kmers.popleft()
			self.bam_records.popleft()
			self.kmer_generators.popleft()
			self.kmer_generators_finished.popleft()
		
	def process_complete_columns(self):
		'''
		Perform calling of columns that are complete, i.e. they cannot receive
		more reads because subsequent reads are further to the right.
		'''
		# compute largest position for which k-mer pileup can safely be generated
		target_pos = self.bam_records[-1].pos + self.k - 1
		self.advance_to(target_pos)
		while self.ref_pos < target_pos:
			yield self.pop_column()


def run_call(reference, bam):
	fasta = pyfaidx.Fasta(reference, as_raw=True)

	samfile = pysam.Samfile(bam)

	k = 8
	window_size = 10
	min_support = 4

	chromosome = None
	for i, bam_alignment in enumerate(samfile):
		if bam_alignment.reference_name != chromosome:
			chromosome = bam_alignment.reference_name
			ref = fasta[chromosome]
#			for h, pos in enumerate_reference_kmers(str(ref).encode('UTF-8'), k):
#				if pos % 10000 == 0:
#					print(h, pos)
#				if pos>10557290:
#					break
			caller = Caller(ref, k)

		caller.add_read(bam_alignment)
		i = 0
		for ref_pos, ref_kmer, pileup_kmers in caller.process_complete_columns():
			if ref_pos % 10000 == 0:
				print(ref_pos)

			#if pileup_kmers != None and (len(pileup_kmers) > 0):
				#print('-'*(k+20))
				#print(hash_to_dna(ref_kmer,k), ref_pos)
				#l = sorted(((count,kmer) for kmer,count in pileup_kmers.items()), reverse=True)
				#for count, kmer in l:
					#print(hash_to_dna(kmer,k), count)
				#i += 1
				##if i == 20:
					##return

			if pileup_kmers != None and (len(pileup_kmers) > 0):
				for kmer,count in pileup_kmers.items():
					if kmer != ref_kmer and count >= min_support:
						print(ref_pos, hash_to_dna(ref_kmer,k), '-->', hash_to_dna(kmer,k), count)

			


		#print(dir(r))
		#print(bam_alignment.cigartuples)
		#print(bam_alignment.get_reference_positions())
		##print(bam_alignment.query_name)
		#print(bam_alignment.query)
		#for i, (h, pos) in enumerate(enumerate_kmers(bam_alignment, k)):
			#print(' '*i, hash_to_dna(h, k), ' ', pos, sep='')
			#print(' '*i, ref[pos-k:pos], ' REF', sep='')
			##pass

		#break
		#if i == 10:
			#break


def main(args):
	run_call(**vars(args))
