"""
Generate sequencing technology specific error profiles
"""
import logging
import pysam
import pyfaidx
from whatshap.core import Caller
from pysam import VariantFile


logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--reference', metavar='FASTA', help='Reference genome', required=True)
    arg('--bam', metavar='BAM', help='Aligned reads', required=True)
    arg('--vcf', metavar='VCF', help='Variants', required=True)
    arg('--kmer', metavar='KMER', help='kmer size', required=True)
    arg('--window', metavar='WINDOW', help='Ignore this many bases on the left and right of each variant position', required=True)


def learn(reference, bam, vcf, kmer, window):
    fasta = pyfaidx.Fasta(reference, as_raw=True)
    bamfile = pysam.AlignmentFile(bam)
    variantslist = []
    call = 0
    vcf_in = VariantFile(vcf)
    for variant in vcf_in.fetch():
        variantslist.append((variant.pos, len(variant.ref)))
    variant = 0
    encoded_references = {}
    chromosome = None
    for bam_alignment in bamfile:
        if not bam_alignment.is_unmapped and not bam_alignment.query_alignment_sequence==None:
            if bam_alignment.reference_name != chromosome:
                chromosome = bam_alignment.reference_name
                if chromosome in encoded_references:
                    caller = Caller(encoded_references[chromosome], int(kmer), int(window))
                else:
                    ref = fasta[chromosome]
                    encoded_references[chromosome] = str(ref).encode('UTF-8')
                    caller = Caller(encoded_references[chromosome], int(kmer), int(window))
            if call==0:
                caller.all_variants(variantslist)
                call = 1
            else:
                pass
            caller.add_read(bam_alignment.pos, bam_alignment.cigartuples, str(bam_alignment.query_alignment_sequence).encode('UTF-8'))
    caller.final_pop()
         


def main(args):
    learn(**vars(args))
	
