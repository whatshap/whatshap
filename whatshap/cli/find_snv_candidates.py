"""
Generate candidate SNP positions.
"""

import pysam
import sys
import re
import pyfaidx
from collections import defaultdict
import datetime
import logging

logger = logging.getLogger(__name__)


# fmt: off
def add_arguments(parser):
    add = parser.add_argument
    add('ref', metavar='REF', help='FASTA with reference genome')
    add('bam', metavar='BAM', help='BAM file')
    add('--minabs', metavar='MIN_ABS', default=3, type=int,
        help='Minimum absolute ALT depth to call a SNP (default: %(default)s).')
    add('--minrel', metavar='MIN_REL', default=0.25, type=float,
        help='Minimum relative ALT depth to call a SNP (default: %(default)s).')
    add('--multi-allelics', default=False, action='store_true',
        help='Also output multi-allelic sites, if not given only the best ALT allele is reported (if unique).')
    add('--sample', metavar='SAMPLE', default='sample',
        help='Put this sample column into VCF (default: output sites-only VCF).')
    add('--chromosome', dest='chromosome', metavar='CHROMOSOME', default=None,
        help='Name of chromosome to process. If not given, all chromosomes are processed.')
    add('-o', '--output', default=sys.stdout, help='Output VCF file.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '--pacbio', dest='datatype', action='store_const', const='pacbio',
        help='Input is PacBio. Sets minrel=0.25 and minabs=3.')
    group.add_argument(
        '--nanopore', dest='datatype', action='store_const', const='nanopore',
        help='Input is Nanopore. Sets minrel=0.4 and minabs=3.')
    group.add_argument(
        '--illumina', dest='datatype', action='store_const', const='illumina',
        help='Input is Illumina. Sets minrel=0.25 and minabs=3.')
# fmt: on


def validate(args, parser):
    pass


def run_find_snv_candidates(
    ref,
    bam,
    minabs=3,
    minrel=0.25,
    multi_allelics=False,
    datatype=None,
    sample="sample",
    chromosome=None,
    output=sys.stdout,
):
    outfile = output
    if output != sys.stdout:
        outfile = open(output, "w")
    if datatype == "pacbio":
        minabs = 3
        minrel = 0.25
    if datatype == "nanopore":
        minabs = 3
        minrel = 0.4
    if datatype == "illumina":
        minabs = 3
        minrel = 0.25
    print(minabs, minrel)
    fasta = pyfaidx.Fasta(ref, as_raw=True)
    print("##fileformat=VCFv4.2", file=outfile)
    print("##fileDate={}".format(datetime.datetime.now().strftime("%Y%m%d")), file=outfile)
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=outfile)
    print('##FILTER=<ID=PASS,Description="All filters passed">', file=outfile)
    header_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if sample is not None:
        header_columns += ["FORMAT", sample]
    print(*header_columns, sep="\t", file=outfile)

    re_nucleotide = re.compile("[ACGTNacgtn]")
    re_indel = re.compile("[-\\+]([0-9]+)")
    re_ref = re.compile("[,\\.]")
    re_ignore = re.compile("([\\$\\*]|\\^.)")

    bamfile = pysam.AlignmentFile(bam, "rb")
    for pileupcolumn in bamfile.pileup(
        contig=chromosome, min_mapping_quality=20, min_base_quality=5
    ):
        try:
            pileup = "".join(
                pileupcolumn.get_query_sequences(mark_matches=True, mark_ends=True, add_indels=True)
            )
        except AssertionError:
            # TODO: seems to happen when there are large insertions/deletions
            continue
        position = pileupcolumn.reference_pos + 1
        chromosome = pileupcolumn.reference_name
        i = 0
        bases = defaultdict(int)
        n = 0
        ref = fasta[chromosome][position - 1].upper()
        if ref == "N":
            continue
        while i < len(pileup):
            m = re_nucleotide.match(pileup[i:])
            if m is not None:
                bases[pileup[i].upper()] += 1
                n += 1
                i += 1
                continue
            m = re_indel.match(pileup[i:])
            if m is not None:
                l = int(m.group(1))
                skip = (m.end() - m.start()) + l
                i += skip
                continue
            m = re_ref.match(pileup[i:])
            if m is not None:
                bases[ref] += 1
                n += 1
                i += 1
                continue
            m = re_ignore.match(pileup[i:])
            if m is not None:
                skip = m.end() - m.start()
                i += skip
                continue
            assert False
        ref_count = bases[ref]
        alts = []
        for base, count in bases.items():
            if base == ref:
                continue
            if (count >= minabs) and (count / (count + ref_count) >= minrel):
                alts.append((count, base))
        alts.sort(reverse=True)
        if len(alts) > 0:
            columns = [chromosome, position, ".", ref, ".", ".", "PASS", "."]
            if sample is not None:
                columns += ["GT", "."]
            if multi_allelics:
                columns[4] = ",".join(base for count, base in alts)
            else:
                # Do we have two equally supported ALT alleles
                if len(alts) > 1 and (alts[0][0] == alts[1][0]):
                    columns[4] = "N"
                    continue
                else:
                    columns[4] = alts[0][1]
            print(*columns, sep="\t", file=outfile)
    if output != sys.stdout:
        outfile.close()


def main(args):
    run_find_snv_candidates(**vars(args))
