"""
Remove phasing information from a VCF file

This script removes all types of phasing information from the input VCF and
prints out the modified VCF to standard output. The modifications are:

- The HP, PS and PQ tags are removed
- Phasing in the GT tag (using pipe notation) is removed. The genotypes are
  sorted in ascending order. For example, a GT value of '1|0' is converted
  to '0/1'.

It is not an error if no phasing information was found.
"""
import sys
import logging
from pysam import VariantFile


logger = logging.getLogger(__name__)

TAGS_TO_REMOVE = frozenset(("HP", "PQ", "PS"))


def add_arguments(parser):
    add = parser.add_argument
    add("vcf", metavar="VCF", help='VCF file. Use "-" to read from standard input')


def unphase_header(header):
    for hr in header.records:
        if hr.key == "phasing":
            hr.remove()
            break

    for tag in TAGS_TO_REMOVE:
        if tag in header.formats:
            header.formats.remove_header(tag)


def run_unphase(vcf_path, outfile):
    """
    Read a VCF file, remove phasing information, and write the result to
    outfile, which must be a file-like object.
    """
    if vcf_path == "-":
        reader = VariantFile(sys.stdin)
    else:
        reader = VariantFile(vcf_path)

    unphase_header(reader.header)
    with VariantFile(outfile, mode="w", header=reader.header) as writer:
        for record in reader:
            for tag in TAGS_TO_REMOVE:
                if tag in record.format:
                    del record.format[tag]
            for call in record.samples.values():
                if (
                    call["GT"] is not None
                    and call["GT"][0] is not None
                    and call["GT"][1] is not None
                ):
                    call["GT"] = sorted(call["GT"])
                call.phased = False
            writer.write(record)


def main(args):
    run_unphase(args.vcf, sys.stdout)
