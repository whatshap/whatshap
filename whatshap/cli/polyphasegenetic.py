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
from whatshap.core import ClusterEditingSolver, Read, ReadSet
from whatshap.cli import log_memory_usage, CommandLineError
from whatshap.polyphaseplots import (
    draw_genetic_clustering,
    draw_genetic_clustering_arrangement,
    draw_phase_comparison,
)
from whatshap.offspringscoring import (
    get_variant_scoring,
    get_phasable_parent_variants,
    get_offspring_gl,
    add_corrected_variant_types,
)
from whatshap.timer import StageTimer
from whatshap.vcf import VcfReader, PhasedVcfWriter, PloidyError
from whatshap.clusterarrangement import arrange_clusters

__author__ = "Sven Schrinner"

PhasingParameter = namedtuple(
    "PhasingParameter",
    ["ploidy", "scoring_window", "simplex", "allow_homozyguous", "allow_deletions"],
)

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
    progeny_file=None,
    ground_truth_file=None,
    scoring_window=160,
    simplex=False,
    allow_homozyguous=False,
    distrust_parent_genotypes=False,
    output=sys.stdout,
    samples=None,
    chromosomes=None,
    indels=True,
    mapping_quality=20,
    tag="PS",
    write_command_line_header=True,
    read_list_filename=None,
    plot=False,
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
        if progeny_file:
            progeny_reader = stack.enter_context(
                VcfReader(
                    progeny_file,
                    indels=indels,
                    phases=True,
                    genotype_likelihoods=False,
                    ploidy=ploidy,
                    mav=True,
                    allele_depth=True,
                )
            )
        if ground_truth_file:
            ground_truth_reader = stack.enter_context(
                VcfReader(
                    ground_truth_file,
                    indels=indels,
                    phases=True,
                    genotype_likelihoods=False,
                    ploidy=ploidy,
                    mav=True,
                    allele_depth=False,
                )
            )

        # exit on non-tetraploid samples
        if ploidy != 4:
            raise CommandLineError(
                "Only ploidy 4 is supported at the moment. Detected was {}".format(ploidy)
            )

        # determine pedigree
        parents, co_parent, offspring = determine_pedigree(
            pedigree_file,
            samples,
            vcf_reader.samples,
            progeny_reader.samples if progeny_file else vcf_reader.samples,
        )

        # validate samples
        vcf_sample_set = set(vcf_reader.samples)
        for sample in samples:
            if sample not in vcf_sample_set:
                raise CommandLineError(
                    "Sample {!r} requested on command-line not found in VCF".format(sample)
                )
        samples = frozenset(samples)

        # Store phasing parameters in tuple to keep function signatures cleaner
        phasing_param = PhasingParameter(
            ploidy=ploidy,
            scoring_window=scoring_window,
            simplex=simplex,
            allow_homozyguous=allow_homozyguous,
            allow_deletions=indels,
        )

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

                logger.info("Number of variants among all samples: %d", len(variant_table))

                # compute scoring matrices for parent samples
                for sample in samples:
                    logger.info("---- Processing individual %s", sample)

                    # compute phasable parent variants
                    varinfo, phasable_indices = get_phasable_parent_variants(
                        variant_table, sample, co_parent[sample], phasing_param
                    )

                    # if progeny file provided, extract region, else just reuse main table
                    timers.start("parse_vcf")
                    if progeny_file:
                        main_positions = [
                            variant_table.variants[i].position for i in phasable_indices
                        ]
                        regions = [
                            (main_positions[i], main_positions[i] + 1)
                            for i in range(len(main_positions))
                        ]
                        progeny_table = progeny_reader.fetch_regions(chromosome, regions)
                    else:
                        progeny_table = variant_table
                    timers.stop("parse_vcf")

                    # compute offspring genotype likelihoods
                    timers.start("scoring")
                    logger.info("Computing genotype likelihoods for offspring ...")
                    if distrust_parent_genotypes:
                        varinfo = add_corrected_variant_types(
                            variant_table,
                            progeny_table,
                            offspring[(sample, co_parent[sample])],
                            varinfo,
                            phasable_indices,
                            phasing_param,
                            0.06,
                        )
                    off_gl, node_to_variant, type_of_node = get_offspring_gl(
                        variant_table,
                        progeny_table,
                        offspring[(sample, co_parent[sample])],
                        varinfo,
                        phasable_indices,
                        phasing_param,
                        0.06,
                    )

                    # store progeny coverage
                    progeny_coverage_all = [0 for _ in range(len(variant_table))]
                    for off in offspring[(sample, co_parent[sample])]:
                        parent_pos = 0
                        progeny_pos = 0
                        allele_depths = progeny_table.allele_depths_of(off)
                        assert len(allele_depths) == len(progeny_table)
                        while progeny_pos < len(allele_depths) and parent_pos < len(variant_table):
                            if (
                                variant_table.variants[parent_pos].position
                                == progeny_table.variants[progeny_pos].position
                            ):
                                progeny_coverage_all[parent_pos] += sum(allele_depths[progeny_pos])
                                progeny_pos += 1
                            else:
                                assert (
                                    variant_table.variants[parent_pos].position
                                    < progeny_table.variants[progeny_pos].position
                                )
                            parent_pos += 1

                    # delete progeny table if dedicated
                    if progeny_file:
                        del progeny_table

                    # compute scoring for variant pairs
                    logger.info("Compute scores for markers ...")
                    scoring = get_variant_scoring(varinfo, off_gl, node_to_variant, phasing_param)

                    # delete genotype likelihoods
                    del off_gl

                    timers.stop("scoring")

                    # cluster variants based on scores
                    timers.start("clustering")
                    logger.info("Clustering marker variants ...")
                    solver = ClusterEditingSolver(scoring, False)
                    clustering = solver.run()
                    del solver

                    # validate clustering
                    for clust in clustering:
                        positions = [node_to_variant[x] for x in clust]
                        if len(set(positions)) != len(clust):
                            print(positions)
                        assert len(set(positions)) == len(clust)
                    timers.stop("clustering")

                    # arrange clusters to haplotypes
                    timers.start("arrangement")
                    logger.info("Arranging clusters ...")
                    haplo_skeletons = arrange_clusters(
                        clustering, node_to_variant, (phasing_param.scoring_window + 1) // 2, ploidy
                    )

                    # determine haplotypes
                    accessible_positions = sorted([v.position for v in variant_table.variants])

                    # for information:
                    # accessible_position maps position (index) of variant_table to genome position
                    # node_to_variant maps a node id of the clustering to a position (index) of the variant_table
                    # varinfo contains the ref- and alt-allele for a position (index) of the variant_table
                    # haplo_skeletons contains ploidy many lists of nodes, which belong to the same haplotype

                    components[sample] = {}
                    superreads[sample] = ReadSet()
                    for i in range(phasing_param.ploidy):
                        superreads[sample].add(Read("superread {}".format(i + 1), 0, 0))

                    signals_per_pos = defaultdict(list)
                    for i, hap in enumerate(haplo_skeletons):
                        for clust in hap:
                            for node in clustering[clust]:
                                signals_per_pos[node_to_variant[node]].append(i)

                    phased_positions = []
                    haplotypes = [[] for _ in range(phasing_param.ploidy)]
                    sample_coverage = []
                    progeny_coverage = []
                    allele_depths = variant_table.allele_depths_of(sample)

                    for pos in range(len(variant_table)):
                        if not phasing_param.allow_homozyguous and len(signals_per_pos[pos]) == 0:
                            continue
                        for i in range(phasing_param.ploidy):
                            if i in signals_per_pos[pos]:
                                allele = varinfo[pos].alt
                            else:
                                allele = varinfo[pos].ref
                            superreads[sample][i].add_variant(accessible_positions[pos], allele, 0)
                            components[sample][accessible_positions[pos]] = accessible_positions[0]
                            haplotypes[i].append(allele)
                        phased_positions.append(accessible_positions[pos])
                        sample_coverage.append(sum(allele_depths[pos]))
                        progeny_coverage.append(progeny_coverage_all[pos])

                    timers.stop("arrangement")

                    timers.start("parse_vcf")
                    if ground_truth_file:
                        regions = [
                            (phased_positions[i], phased_positions[i] + 1)
                            for i in range(len(phased_positions))
                        ]
                        ground_truth_table = ground_truth_reader.fetch_regions(chromosome, regions)
                    timers.stop("parse_vcf")

                    # create plots
                    if plot:
                        num_vars = max([max(c) for c in clustering])
                        draw_genetic_clustering(clustering, num_vars, output + ".clusters.pdf")

                        draw_genetic_clustering_arrangement(
                            clustering,
                            node_to_variant,
                            haplo_skeletons,
                            type_of_node,
                            (phasing_param.scoring_window + 1) // 2,
                            num_vars,
                            output + ".arrangement.pdf",
                        )

                        if ground_truth_file:
                            draw_phase_comparison(
                                haplotypes,
                                phased_positions,
                                progeny_coverage,
                                variant_table,
                                ground_truth_table,
                                output + ".comparison.pdf",
                            )

                with timers("write_vcf"):
                    logger.info("======== Writing VCF")
                    vcf_writer.write(
                        chromosome,
                        superreads,
                        components,
                    )
                    logger.info("Done writing VCF")
                logger.debug("Chromosome %r finished", chromosome)
                timers.start("parse_vcf")
            timers.stop("parse_vcf")
        except PloidyError as e:
            raise CommandLineError(e)

    logger.info("\n== SUMMARY ==")

    log_memory_usage()
    logger.info("Time spent parsing VCF:                   %6.1f s", timers.elapsed("parse_vcf"))
    logger.info("Time spent for genetic scoring:           %6.1f s", timers.elapsed("scoring"))
    logger.info("Time spent for clustering:                %6.1f s", timers.elapsed("clustering"))
    logger.info("Time spent for cluster arrangement:       %6.1f s", timers.elapsed("arrangement"))
    logger.info("Time spent writing VCF:                   %6.1f s", timers.elapsed("write_vcf"))
    logger.info("Time spent on rest:                       %6.1f s", timers.total() - timers.sum())
    logger.info("Total elapsed time:                       %6.1f s", timers.total())


def determine_pedigree(pedigree_file, samples, vcf_samples, progeny_samples):

    parents, co_parent, offspring = dict(), dict(), defaultdict(list)
    with open(pedigree_file, "r") as ped:
        for line in ped:
            tokens = line.replace("\n", "").split(" ")
            if len(tokens) != 3:
                logger.error("Malformed pedigree file: {}".format(line))
                raise CommandLineError(None)
            for token in tokens[:2]:
                if token not in vcf_samples:
                    logger.error(
                        "Parent sample {} from pedigree file is not present in main VCF file".format(
                            token
                        )
                    )
                    raise CommandLineError(None)
            if tokens[2] not in progeny_samples:
                logger.error(
                    "Progeny sample {} from pedigree file is not present in progeny (or main) VCF file".format(
                        token[2]
                    )
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
                    "Sample {} does not have a co-parent for the pedigree phasing".format(sample)
                )
                raise CommandLineError(None)
            if len(offspring[(sample, co_parent[sample])]) == 0:
                logger.error(
                    "Sample {} does not have any offspring according to pedigree file".format(
                        sample
                    )
                )
                raise CommandLineError(None)

    return parents, co_parent, offspring


def add_arguments(parser):
    arg = parser.add_argument
    # Positional argument
    arg(
        "variant_file",
        metavar="VCF",
        help="VCF file with variants to be phased (can be gzip-compressed)",
    )
    arg("pedigree_file", metavar="PEDIGREE", help="Pedigree file.")
    arg(
        "-P",
        "--progeny_file",
        required=False,
        help="File with progeny genotypes. If not specified, information is taken from main input file.",
    )
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
    arg(
        "--ground-truth-file",
        "-g",
        required=False,
        help="File containing a ground truth phasing. Only used for plotting.",
    )
    arg(
        "--scoring-window",
        metavar="SCORINGWINDOW",
        dest="scoring_window",
        type=int,
        default=160,
        required=False,
        help="Size of the window (in variants) for statistical offspring scoring.",
    )
    arg(
        "--simplex",
        dest="simplex",
        default=False,
        action="store_true",
        help="Reduce offspring scoring to simplex-nulliplex variants only.",
    )
    arg(
        "--allow-homozyguous",
        dest="allow_homozyguous",
        default=False,
        action="store_true",
        help="Writes sides which are phased as homozyguous into output instead of old genotype.",
    )
    arg(
        "--distrust-parent-genotypes",
        dest="distrust_parent_genotypes",
        default=False,
        action="store_true",
        help="Uses progeny genotypes to double-check parental genotypes.",
    )
    arg(
        "--plot",
        dest="plot",
        default=False,
        action="store_true",
        help="Plots the variant clustering and arrangement.",
    )


def validate(args, parser):
    pass


def main(args):
    run_polyphasegenetic(**vars(args))
