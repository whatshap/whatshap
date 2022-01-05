"""
Call variants
"""
import logging
import pysam
import pyfaidx
from collections import defaultdict, namedtuple, deque
from contextlib import ExitStack
from whatshap.core import Caller
from pysam import VariantFile
#from whatshap._call import hash_to_dna
#from whatshap.align import edit_distance
#import cProfile
#import re
#import timeit

logger = logging.getLogger(__name__)


def add_arguments(parser):
	add = parser.add_argument
	add('reference', metavar='FASTA', help='Reference genome FASTA')
	add('bam', metavar='BAM', help='BAM file ')
	add('vcf', metavar='VCF', help='VCF file ')
	add('kmer_size', metavar='KMER', help='kmer_size')
	add('window_size', metavar='WINDOW', help='length of the window in one direction i.e. half the total window size')
	add('epsilon', metavar='EPSILON', help='value used for combination we never saw')


def run_call(reference, bam, vcf, kmer_size, window_size, epsilon):
    fasta = pyfaidx.Fasta(reference, as_raw=True)
    bamfile = pysam.AlignmentFile(bam)
    variantslist=[]
    #reader= open(vcf, 'r')
    call=0
    vcf_in = VariantFile(vcf)
    for variant in vcf_in.fetch():
        #line_r = line.strip().split("\t")
        #varpos= int(line_r[1])
        variantslist.append((variant.pos, len(variant.ref)))
    variant=0
    encoded_references={}
    k = int(kmer_size)
    window= int(window_size)
    epsilon= float(epsilon)
    chromosome = None
    for bam_alignment in bamfile:
        if not bam_alignment.is_unmapped and not bam_alignment.query_alignment_sequence==None:
            if bam_alignment.reference_name != chromosome:
            ###here we can add a line or two checking if this is not the first chromosome then the remaining kmers form the previous reference need to be flushed out
            #depends on if final_pop function works as expected
                chromosome = bam_alignment.reference_name
                if chromosome in encoded_references:
                    caller = Caller(encoded_references[chromosome], k, epsilon, window)
                else:
                    ref = fasta[chromosome]
                    #caller = Caller(str(ref).encode('UTF-8'), k)
                    encoded_references[chromosome]= str(ref).encode('UTF-8')

                    caller = Caller(encoded_references[chromosome], k, epsilon, window)
            if call==0:
                caller.all_variants(variantslist)
                call=1
            else:
                pass
            caller.add_read(bam_alignment.pos, bam_alignment.cigartuples, str(bam_alignment.query_alignment_sequence).encode('UTF-8'))
    caller.final_pop()
    caller.calc_probs()
         


def main(args):
	#cProfile.runctx('run_call(**vars(args))', {'run_call': run_call, 'args':args}, {})
	#code_to_test = """
	run_call(**vars(args))
	#"""
	#global argms
	#argms= args
	#elapsed_time = timeit.timeit(stmt= "run_call(**vars(argms))",globals=globals(), number=100)/100
	#print(elapsed_time)
