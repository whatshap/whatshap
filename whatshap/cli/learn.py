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
    arg("bam", metavar="BAM", help="Read alignments")
    arg("vcf", metavar="VCF", help="List of variants")
    arg("--reference", "-r", metavar="FASTA", help="Reference genome", required=True)
    arg("-k", "--kmer", dest="k", metavar="K", help="k-mer size", type=int, default=7)
    arg(
        "--window",
        "-w",
        metavar="WINDOW",
        help="Ignore this many bases on the left and right of each variant position",
        type=int,
        default=25,
    )
    arg("--output", "-o", metavar="OUT", help="Output file with kmer-pair counts", required=True)


def run_learn(reference, bam, vcf, k: int, window: int, output):
    with VariantFile(vcf) as vcf:
        variants = [(variant.pos, len(variant.ref)) for variant in vcf.fetch()]

    with pyfaidx.Fasta(reference, as_raw=True) as fasta, pysam.AlignmentFile(bam) as bamfile:
        call = 0
        encoded_references = {}
        chromosome = None
        open(output, "w").close()
        output_c = str(output).encode("UTF-8")
        for alignment in bamfile:
            if alignment.is_unmapped or alignment.query_alignment_sequence is None:
                continue
            if alignment.reference_name != chromosome:
                chromosome = alignment.reference_name
                if chromosome not in encoded_references:
                    ref = fasta[chromosome]
                    encoded_references[chromosome] = str(ref).encode("UTF-8")
                caller = Caller(encoded_references[chromosome], k, window)
            if call == 0:
                caller.all_variants(variants)
                call = 1
            caller.add_read(
                alignment.pos,
                alignment.cigartuples,
                str(alignment.query_alignment_sequence).encode("UTF-8"),
                output_c,
            )
        caller.final_pop(output_c)


def main(args):
    run_learn(**vars(args))
