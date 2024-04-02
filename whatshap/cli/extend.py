"""
Extend phasing information from haplotagged reads to variants
"""

import itertools
import logging
import sys
from contextlib import ExitStack
from typing import List, Optional, Union, Dict, Tuple


from whatshap import __version__
from whatshap.cli import PhasedInputReader, CommandLineError, log_memory_usage
from whatshap.core import NumericSampleIds, Variant, Read
from whatshap.timer import StageTimer
from whatshap.utils import IndexedFasta
from whatshap.vcf import VcfReader, PhasedVcfWriter, VcfError

logger = logging.getLogger(__name__)


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    arg("-o", "--output",
        default=sys.stdout,
        help="Output file. If omitted, use standard output.")
    arg("--reference", "-r", metavar="FASTA",
        help="Reference file. Must be accompanied by .fai index (create with samtools faidx)")
    arg("--gap-threshold", "-g", metavar="GAPTHRESHOLD", default=70, type=int,
        help="Threshold percentage of qualities for assigning phase information to a variant.")
    arg("--cut-poly", "-c", metavar="CUTPOLY", default=10, type=int,
        help="ignore polymers longer than the cut value.")
    arg("--only-indels", "-i", default=False, action="store_true",
        help="extend new phasing information only to indels.")
    arg("--ignore-read-groups", default=False, action="store_true",
        help="Ignore read groups in BAM/CRAM header and assume all reads come from the same sample.")
    arg("--chromosome", dest="chromosomes", metavar="CHROMOSOME", default=[], action="append",
        help="Name of chromosome to phase. If not given, all chromosomes in the input VCF are phased. "
        "Can be used multiple times.")
    arg("variant_file", metavar="VCF", help="VCF file with phased variants (must be gzip-compressed and indexed)")
    arg("alignment_file", metavar="ALIGNMENTS",
        help="BAM/CRAM file with alignments to be tagged by haplotype")


# fmt: on


def run_extend(
    variant_file,
    alignment_file,
    output=None,
    reference: Union[None, bool, str] = False,
    ignore_read_groups: bool = False,
    only_indels: bool = False,
    chromosomes: Optional[List[str]] = None,
    gap_threshold: int = 70,
    cut_poly: int = 10,
    write_command_line_header: bool = True,
    tag: str = "PS",
):
    if reference is None:
        raise CommandLineError("Option --reference should be specified")
    timers = StageTimer()
    timers.start("extend-run")
    command_line: Optional[str]
    if write_command_line_header:
        command_line = "(whatshap {}) {}".format(__version__, " ".join(sys.argv[1:]))
    else:
        command_line = None
    with ExitStack() as stack:
        logger.debug("Creating PhasedInputReader")
        phased_input_reader = stack.enter_context(
            PhasedInputReader(
                [alignment_file],
                None if reference is False else reference,
                NumericSampleIds(),
                ignore_read_groups,
                only_snvs=False,
            )
        )
        logger.debug("Creating PhasedVcfWriter")
        try:
            vcf_writer = stack.enter_context(
                PhasedVcfWriter(
                    command_line=command_line,
                    in_path=variant_file,
                    out_file=output,
                    tag=tag,
                )
            )
        except (OSError, VcfError) as e:
            raise CommandLineError(e)

        vcf_reader = stack.enter_context(VcfReader(variant_file, phases=True))

        if ignore_read_groups and len(vcf_reader.samples) > 1:
            raise CommandLineError(
                "When using --ignore-read-groups on a VCF with "
                "multiple samples, --sample must also be used."
            )
        with timers("read-fasta"):
            fasta = stack.enter_context(IndexedFasta(reference))
        for variant_table in timers.iterate("parse-vcf", vcf_reader):
            chromosome = variant_table.chromosome
            fasta_chr = fasta[chromosome]
            logger.info(f"Processing chromosome {chromosome}...")
            if chromosomes and chromosome not in chromosomes:
                logger.info(
                    f"Leaving chromosome {chromosome} unchanged "
                    "(present in VCF, but not requested by --chromosome)"
                )
                with timers("write-vcf"):
                    vcf_writer.write_unchanged(chromosome)
                continue
            sample_to_super_reads, sample_to_components = (dict(), dict())
            for sample in vcf_reader.samples:
                logger.info(f"Processing sample {sample}")
                with timers("read-bam"):
                    reads, _ = phased_input_reader.read(chromosome, variant_table.variants, sample)
                phases = variant_table.phases_of(sample)
                genotypes = variant_table.genotypes_of(sample)
                homozygous = dict()
                change = dict()
                phased = dict()
                homozygous_number = 0
                phased_number = 0
                for variant, (phase, genotype) in zip(
                    variant_table.variants, zip(phases, genotypes)
                ):
                    homozygous[variant.position] = genotype.is_homozygous()
                    phased[variant.position] = phase
                    phased_number += phase is not None
                    homozygous_number += genotype.is_homozygous()
                    change[variant.position] = variant
                logger.info(f"Number of homozygous variants is {homozygous_number}")
                logger.info(f"Number of already phased variants is {phased_number}")
                with timers("compute-votes"):
                    votes = compute_votes(homozygous, reads)
                with timers("compute-consensus"):
                    sample_to_super_reads[sample], sample_to_components[sample] = consensus(
                        only_indels, gap_threshold, cut_poly, fasta_chr, change, phased, votes
                    )
            with timers("write-vcf"):
                vcf_writer.write(chromosome, sample_to_super_reads, sample_to_components)
    timers.stop("extend-run")
    log_time_and_memory_usage(timers)


def log_time_and_memory_usage(timers):
    logger.info("\n# Resource usage")
    log_memory_usage()
    # fmt: off
    logger.info("Finished in :                              %6.1f s", timers.elapsed("extend-run"))
    logger.info("Time spent reading reference:              %6.1f s", timers.elapsed("read-fasta"))
    logger.info("Time spent reading VCF:                    %6.1f s", timers.elapsed("parse-vcf"))
    logger.info("Time spent writing VCF:                    %6.1f s", timers.elapsed("write-vcf"))
    logger.info("Time spent reading BAM:                    %6.1f s", timers.elapsed("read-bam"))
    logger.info("Time spent computing votes:                %6.1f s", timers.elapsed("compute-votes"))
    logger.info("Time spent spent computing consensus:      %6.1f s", timers.elapsed("compute-consensus"))
    # fmt: on


def consensus(only_indels, gap_threshold, cut_poly, fasta_chr, change, phased, votes):
    super_reads = [[], []]
    components = dict()

    for pos, var in votes.items():
        al1, q, score1 = best_candidate(components, pos, var)
        if phased[pos] is None:
            if 100 * q < gap_threshold:
                continue
            if only_indels and change[pos].is_snv():
                continue
            if cut_poly > 0:
                max_length = max(
                    length_of_polymer(fasta_chr, pos + 1, 1, cut_poly),
                    length_of_polymer(fasta_chr, pos, -1, cut_poly),
                )
                if max_length >= cut_poly:
                    continue
        super_reads[0].append(Variant(pos, allele=al1, quality=score1))
        super_reads[1].append(Variant(pos, allele=al1 ^ 1, quality=score1))
    for read in super_reads:
        read.sort(key=lambda x: x.position)
    return super_reads, components


def best_candidate(components, pos, var):
    lst = list(var.items())
    lst.sort(key=lambda x: x[-1], reverse=True)
    (ps1, al1), score1 = lst[0]
    total = sum(e[-1] for e in lst)
    components[pos] = ps1
    q = score1 / total
    return al1, q, score1


def length_of_polymer(ref: str, start: int, step: int, threshold: int) -> int:
    res = 0
    for i in itertools.count(start, step):
        if res < threshold and ref[i] == ref[start]:
            res += 1
        else:
            break
    return res


def compute_votes(
    homozygous: Dict[int, bool],
    reads: List[Read],
) -> Dict[int, Dict[Tuple[int, int], int]]:
    votes = dict()
    number_of_skipped = 0
    for read in reads:
        ps, ht = read.PS_tag - 1, read.HP_tag - 1
        if ht < 0 or ps < 0:
            continue
        if ht > 1:
            number_of_skipped += 1
        for variant in read:
            if homozygous[variant.position]:
                continue
            if variant.position not in votes:
                votes[variant.position] = dict()
            if (ps, 0) not in votes[variant.position]:
                votes[variant.position][(ps, 0)] = 0
                votes[variant.position][(ps, 1)] = 0
            votes[variant.position][(ps, ht ^ variant.allele)] += variant.quality
    if number_of_skipped > 0:
        logger.warning(
            f"{number_of_skipped} reads were skipped due incorrect HP. The extend command supports only a diploid input"
        )
    return votes


def main(args):
    run_extend(**vars(args))
