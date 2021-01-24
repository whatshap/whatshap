"""
Phase variants in a polyploid VCF using a clustering+threading algorithm.

Read a VCF and one or more files with phase information (BAM/CRAM or VCF phased
blocks) and phase the variants. The phased VCF is written to standard output.
Requires to specify a ploidy for the phasable input. Allows to specify a block
cut sensitivity to balance out length and accuracy of phased blocks.

"""
import sys
import logging
import platform

from collections import namedtuple, defaultdict

from contextlib import ExitStack

from whatshap import __version__
from whatshap.core import ClusterEditingSolver
from whatshap.cli import log_memory_usage, CommandLineError
from whatshap.polyphaseplots import draw_genetic_clustering
from whatshap.offspringscoring import get_variant_scoring
from whatshap.timer import StageTimer
from whatshap.vcf import VcfReader, PhasedVcfWriter, PloidyError

__author__ = "Sven Schrinner"

PhasingParameter = namedtuple("PhasingParameter", ["ploidy",],)

logger = logging.getLogger(__name__)


def print_readset(readset):
    result = ""
    positions = readset.get_positions()
    for read in readset:
        result += read.name + "\t" + "\t" + "\t"
        for pos in positions:
            if pos in read:
                # get corresponding variant
                for var in read:
                    if var.position == pos:
                        result += str(var.allele)
            else:
                result += " "
        result += "\n"
    print(result)


def run_polyphasegenetic(
    variant_file,
    pedigree_file,
    ploidy,
    output=sys.stdout,
    samples=None,
    chromosomes=None,
    indels=True,
    mapping_quality=20,
    tag="PS",
    write_command_line_header=True,
    read_list_filename=None,
):
    """
    Run Polyploid Phasing.

    variant-file -- path to input VCF
    pedigree_file -- pedigree file
    output -- path to output VCF or a file like object
    samples -- names of samples to phase. An empty list means: phase all samples
    chromosomes -- names of chromosomes to phase. An empty list means: phase all chromosomes
    ignore_read_groups
    mapping_quality -- discard reads below this mapping quality
    tag -- How to store phasing info in the VCF, can be 'PS' or 'HP'
    write_command_line_header -- whether to add a ##commandline header to the output VCF
    """
    timers = StageTimer()
    logger.info(
        "This is WhatsHap (polyploid-genetic) %s running under Python %s",
        __version__,
        platform.python_version(),
    )
    with ExitStack() as stack:

        if write_command_line_header:
            command_line = "(whatshap {}) {}".format(__version__, " ".join(sys.argv[1:]))
        else:
            command_line = None
        try:
            vcf_writer = stack.enter_context(
                PhasedVcfWriter(
                    command_line=command_line,
                    in_path=variant_file,
                    out_file=output,
                    tag=tag,
                    ploidy=ploidy,
                )
            )
        except OSError as e:
            raise CommandLineError(e)

        vcf_reader = stack.enter_context(
            VcfReader(
                variant_file,
                indels=indels,
                phases=True,
                genotype_likelihoods=False,
                ploidy=ploidy,
                mav=True,
                allele_depth=True,
            )
        )

        # exit on non-tetraploid samples
        if ploidy != 4:
            raise CommandLineError("Only ploidy 4 is supported. Detected was {}".format(ploidy))

        # determine pedigree
        parents, co_parent, offspring = dict(), dict(), defaultdict(list)

        with open(pedigree_file, "r") as ped:
            for line in ped:
                tokens = line.replace("\n", "").split(" ")
                if len(tokens) != 3:
                    logger.error("Malformed pedigree file: {}".format(line))
                    raise CommandLineError(None)
                for token in tokens:
                    if token not in vcf_reader.samples:
                        logger.error(
                            "Sample {} from pedigree file is not present in VCF file".format(token)
                        )
                        raise CommandLineError(None)

                if tokens[2] in parents:
                    logger.error(
                        "Sample {} from pedigree file is listed as offspring multiple times".format(
                            tokens[2]
                        )
                    )
                    raise CommandLineError(None)
                if tokens[0] in co_parent and co_parent[tokens[0]] != tokens[1]:
                    logger.error(
                        "Sample {} from pedigree file has multiple co-parents".format(tokens[0])
                    )
                    raise CommandLineError(None)
                if tokens[1] in co_parent and co_parent[tokens[1]] != tokens[0]:
                    logger.error(
                        "Sample {} from pedigree file has multiple co-parents".format(tokens[1])
                    )
                    raise CommandLineError(None)

                co_parent[tokens[0]] = tokens[1]
                co_parent[tokens[1]] = tokens[0]
                parents[tokens[2]] = (tokens[0], tokens[1])
                offspring[(tokens[0], tokens[1])].append(tokens[2])
                offspring[(tokens[1], tokens[0])].append(tokens[2])
                # print("trio: {} + {} = {}".format(tokens[0], tokens[1], tokens[2]))

        if not samples:
            samples = [parent for parent in co_parent]
        else:
            for sample in samples:
                if sample not in co_parent:
                    logger.error(
                        "Sample {} does not have a co-parent for the pedigree phasing".format(
                            sample
                        )
                    )
                    raise CommandLineError(None)
                if len(offspring[(sample, co_parent[sample])]) == 0:
                    logger.error(
                        "Sample {} does not have any offspring according to pedigree file".format(
                            sample
                        )
                    )
                    raise CommandLineError(None)

        vcf_sample_set = set(vcf_reader.samples)
        for sample in samples:
            if sample not in vcf_sample_set:
                raise CommandLineError(
                    "Sample {!r} requested on command-line not found in VCF".format(sample)
                )

        samples = frozenset(samples)

        # Store phasing parameters in tuple to keep function signatures cleaner
        phasing_param = PhasingParameter(ploidy=ploidy,)

        timers.start("parse_vcf")
        try:
            for variant_table in vcf_reader:
                chromosome = variant_table.chromosome
                timers.stop("parse_vcf")
                if (not chromosomes) or (chromosome in chromosomes):
                    logger.info("======== Working on chromosome %r", chromosome)
                else:
                    logger.info(
                        "Leaving chromosome %r unchanged (present in VCF but not requested by option --chromosome)",
                        chromosome,
                    )
                    with timers("write_vcf"):
                        superreads, components = dict(), dict()
                        vcf_writer.write(chromosome, superreads, components)
                    continue

                # These two variables hold the phasing results for all samples
                superreads, components = dict(), dict()

                logger.info(
                    "Number of variants among all samples: %d", len(variant_table),
                )

                # compute scoring matrices for parent samples
                for sample in samples:
                    logger.info("---- Processing individual %s", sample)
                    timers.start("scoring")
                    scoring, node_to_variant = get_variant_scoring(
                        variant_table,
                        sample,
                        co_parent[sample],
                        offspring[(sample, co_parent[sample])],
                        phasing_param,
                    )
                    timers.stop("scoring")

                    timers.start("phasing")
                    solver = ClusterEditingSolver(scoring, False)
                    clustering = solver.run()
                    del solver

                    with open(output+".clusters.txt", "w") as out:
                        var_to_position = [var.position+1 for var in variant_table.variants]

                        for i, cluster in enumerate(sorted(clustering, key=lambda x: -len(x))):
                            out.write("Cluster {}: {}".format(i, " ".join(list(map(lambda x: str(var_to_position[node_to_variant[x]]), cluster)))))
                            out.write("\n")
                    
                    timers.stop("phasing")

                    num_vars = max([max(c) for c in clustering])
                    draw_genetic_clustering(
                        clustering, num_vars, output + ".clusters.pdf",
                    )

                with timers("write_vcf"):
                    logger.info("======== Writing VCF")
                    """
                    vcf_writer.write(
                        chromosome,
                        superreads,
                        components,
                    )
                    """
                    # TODO: Use genotype information to polish results
                    # assert len(changed_genotypes) == 0
                    logger.info("Done writing VCF")
                logger.debug("Chromosome %r finished", chromosome)
                timers.start("parse_vcf")
            timers.stop("parse_vcf")
        except PloidyError as e:
            raise CommandLineError(e)

    logger.info("\n== SUMMARY ==")

    log_memory_usage()
    logger.info(
        "Time spent reading BAM/CRAM:                 %6.1f s", timers.elapsed("read_bam"),
    )
    logger.info(
        "Time spent parsing VCF:                      %6.1f s", timers.elapsed("parse_vcf"),
    )
    logger.info(
        "Time spent for genetic scoring:              %6.1f s", timers.elapsed("scoring"),
    )
    logger.info(
        "Time spent for genetic phasing:              %6.1f s", timers.elapsed("phasing"),
    )
    logger.info(
        "Time spent writing VCF:                      %6.1f s", timers.elapsed("write_vcf"),
    )
    logger.info(
        "Time spent on rest:                          %6.1f s", timers.total() - timers.sum(),
    )
    logger.info(
        "Total elapsed time:                          %6.1f s", timers.total(),
    )


def add_arguments(parser):
    arg = parser.add_argument
    # Positional argument
    arg(
        "variant_file",
        metavar="VCF",
        help="VCF file with variants to be phased (can be gzip-compressed)",
    )
    arg(
        "pedigree_file", metavar="PEDIGREE", help="Pedigree file.",
    )
    # arg('ploidy', metavar='PLOIDY', type=int,
    #    help='The ploidy of the sample(s).')

    arg(
        "-o",
        "--output",
        default=sys.stdout,
        help="Output VCF file. Add .gz to the file name to get compressed output. "
        "If omitted, use standard output.",
    )
    arg(
        "--tag",
        choices=("PS", "HP"),
        default="PS",
        help="Store phasing information with PS tag (standardized) or "
        "HP tag (used by GATK ReadBackedPhasing) (default: %(default)s)",
    )

    arg = parser.add_argument_group("Input pre-processing, selection, and filtering").add_argument
    arg(
        "--indels",
        dest="indels",
        default=False,
        action="store_true",
        help="Also phase indels (default: do not phase indels)",
    )
    arg(
        "--sample",
        dest="samples",
        metavar="SAMPLE",
        default=[],
        action="append",
        help="Name of a sample to phase. If not given, all samples in the "
        "input VCF are phased. Can be used multiple times.",
    )
    arg(
        "--chromosome",
        dest="chromosomes",
        metavar="CHROMOSOME",
        default=[],
        action="append",
        help="Name of chromosome to phase. If not given, all chromosomes in the "
        "input VCF are phased. Can be used multiple times.",
    )

    # add polyphase specific arguments
    arg = parser.add_argument_group("Parameters for phasing steps").add_argument
    arg(
        "--ploidy",
        "-p",
        metavar="PLOIDY",
        type=int,
        required=True,
        help="The ploidy of the sample(s). Argument is required.",
    )


def validate(args, parser):
    pass


def main(args):
    run_polyphasegenetic(**vars(args))
