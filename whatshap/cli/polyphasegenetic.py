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
import argparse

from collections import defaultdict
from dataclasses import dataclass

from contextlib import ExitStack

from whatshap import __version__
from whatshap.core import Read, ReadSet
from whatshap.cli import log_memory_usage, CommandLineError
from whatshap.polyphase.clusterarrangement import arrange_clusters
from whatshap.polyphase.plots import create_genetic_plots
from whatshap.polyphase.solver import ClusterEditingSolver
from whatshap.polyphase.variantselection import compute_phasable_variants, filter_variants
from whatshap.polyphase.offspringscoring import (
    get_offspring_gl,
    correct_variant_types,
    get_variant_scoring,
)
from whatshap.timer import StageTimer
from whatshap.vcf import VcfReader, PhasedVcfWriter, PloidyError

__author__ = "Sven Schrinner"


@dataclass
class PolyphaseGeneticParameter:
    ploidy: int
    scoring_window: int
    allele_error_rate: float
    complexity_support: int
    ratio_cutoff: float
    distrust_genotypes: bool
    allow_deletions: bool
    plot: bool
    output: str


logger = logging.getLogger(__name__)


def run_polyphasegenetic(
    variant_file,
    pedigree_file,
    ploidy,
    progeny_file=None,
    ground_truth_file=None,
    scoring_window=250,
    allele_error_rate=0.06,
    ratio_cutoff=0.0,
    complexity_support=0,
    distrust_genotypes=False,
    output=sys.stdout,
    samples=None,
    chromosomes=None,
    only_snvs=False,
    tag="PS",
    write_command_line_header=True,
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
                    mav=False,
                )
            )
        except OSError as e:
            raise CommandLineError(e)

        parent_reader = stack.enter_context(
            VcfReader(
                variant_file,
                only_snvs=only_snvs,
                phases=True,
                genotype_likelihoods=False,
                ploidy=ploidy,
                mav=True,
                allele_depth=True,
            )
        )
        progeny_reader = None
        if progeny_file:
            progeny_reader = stack.enter_context(
                VcfReader(
                    progeny_file,
                    only_snvs=only_snvs,
                    phases=True,
                    genotype_likelihoods=False,
                    ploidy=ploidy,
                    mav=True,
                    allele_depth=True,
                )
            )

        # determine pedigree
        parent_file_samples = parent_reader.samples
        progeny_file_samples = progeny_reader.samples if progeny_reader else None
        samples, sample_to_coparent, sample_to_progeny = determine_pedigree(
            pedigree_file, samples, parent_file_samples, progeny_file_samples
        )

        # validate samples
        parent_sample_set = set(parent_reader.samples)
        for sample in samples:
            if sample not in parent_sample_set:
                raise CommandLineError(
                    "Sample {!r} requested on command-line not found in VCF".format(sample)
                )
        samples = frozenset(samples)

        # Store phasing parameters in tuple to keep function signatures cleaner
        phasing_param = PolyphaseGeneticParameter(
            ploidy=ploidy,
            scoring_window=scoring_window,
            allele_error_rate=allele_error_rate,
            complexity_support=complexity_support,
            ratio_cutoff=ratio_cutoff,
            distrust_genotypes=distrust_genotypes,
            allow_deletions=not only_snvs,
            plot=plot,
            output=output,
        )

        timers.start("parse_vcf")
        try:
            for variant_table in parent_reader:
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

                # Phase one sample at a time
                for sample in samples:
                    logger.info("---- Processing individual %s", sample)
                    coparent = sample_to_coparent[sample]
                    progeny_list = sample_to_progeny[sample]
                    logger.info("Detected %s as co-parent for %s.", coparent, sample)

                    superreads[sample], components[sample] = phase_single_sample(
                        chromosome,
                        progeny_reader,
                        ground_truth_file,
                        sample,
                        coparent,
                        progeny_list,
                        variant_table,
                        timers,
                        phasing_param,
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
    if plot:
        logger.info("Time spent creating plots:                %6.1f s", timers.elapsed("plots"))
    logger.info("Time spent on rest:                       %6.1f s", timers.total() - timers.sum())
    logger.info("Total elapsed time:                       %6.1f s", timers.total())


def phase_single_sample(
    chromosome,
    progeny_reader,
    ground_truth_reader,
    sample,
    coparent,
    progeny_list,
    variant_table,
    timers,
    param,
):
    # compute phasable parent variants
    varinfo = compute_phasable_variants(variant_table, sample, coparent, param)

    # if progeny file provided, extract region, else just reuse main table
    timers.start("parse_vcf")
    logger.info("Extracting progeny allele depths ...")
    if progeny_reader:
        positions = [variant_table.variants[i].position for i in varinfo.get_phasable()]
        regions = [(positions[i], positions[i] + 1) for i in range(len(positions))]
        progeny_table = progeny_reader.fetch_regions(chromosome, regions)
    else:
        progeny_table = variant_table
    timers.stop("parse_vcf")

    # store progeny coverage
    parent_cov, co_parent_cov, progeny_cov = get_parent_progeny_coverage(
        sample, coparent, progeny_list, variant_table, progeny_table
    )

    # filter variants based on coverage ratio between parent, co-parent and progeny
    if param.ratio_cutoff > 1.0:
        logger.info("Filtering variant positions based on coverage ratios ...")
        old_num = len(varinfo.get_phasable())
        filter_variants(varinfo, parent_cov, co_parent_cov, progeny_cov, param.ratio_cutoff)
        logger.info("Kept %i out of %i variants.", len(varinfo.get_phasable()), old_num)

    # compute offspring genotype likelihoods
    timers.start("scoring")
    logger.info("Computing progeny genotype likelihoods ...")
    if param.distrust_genotypes:
        correct_variant_types(variant_table, progeny_table, progeny_list, varinfo, param)
    off_gl = get_offspring_gl(variant_table, progeny_table, progeny_list, varinfo, param)

    # delete progeny table if dedicated to reduce memory footprint
    if progeny_reader:
        del progeny_table

    # compute scoring for variant pairs
    logger.info("Compute scores for markers ...")
    scoring = get_variant_scoring(varinfo, off_gl, param)

    # delete genotype likelihoods to reduce memory footprint
    del off_gl

    timers.stop("scoring")

    # cluster variants based on scores
    timers.start("clustering")
    logger.info("Clustering marker alleles ...")
    solver = ClusterEditingSolver(scoring, False)
    clustering = solver.run()
    del solver
    assert clustering
    assert any(len(c) > 1 for c in clustering)
    timers.stop("clustering")

    # arrange clusters to haplotypes
    timers.start("arrangement")
    logger.info("Arranging clusters ...")
    padding = int(param.scoring_window * 3.0 + 1)
    haplo_skeletons = arrange_clusters(clustering, padding, param.ploidy)

    # determine haplotypes
    accessible_positions = sorted([v.position for v in variant_table.variants])

    # for information:
    # accessible_position maps position (index) of variant_table to genome position
    # varinfo.node_to_variant() maps a node id to a position (index) of the variant_table
    # varinfo contains the ref- and alt-allele for a position (index) of the variant_table
    # haplo_skeletons contains ploidy many lists of nodes, which belong to the same haplotype

    components = {}
    superreads = ReadSet()
    for i in range(param.ploidy):
        superreads.add(Read("superread {}".format(i + 1), 0, 0))

    marker_per_pos = defaultdict(list)
    for i, hap in enumerate(haplo_skeletons):
        for clust in hap:
            for node in clustering[clust]:
                marker_per_pos[varinfo.node_to_variant(node)].append(i)

    phased_positions = []
    haplotypes = [[] for _ in range(param.ploidy)]
    parent_coverage = []
    co_parent_coverage = []
    progeny_coverage = []

    for pos in range(len(variant_table)):
        if len(marker_per_pos[pos]) == 0:
            continue
        for i in range(param.ploidy):
            if i in marker_per_pos[pos]:
                allele = varinfo[pos].alt
            else:
                allele = varinfo[pos].ref
            superreads[i].add_variant(accessible_positions[pos], allele, 0)
            components[accessible_positions[pos]] = accessible_positions[0]
            haplotypes[i].append(allele)
        phased_positions.append(accessible_positions[pos])
        parent_coverage.append(parent_cov[pos])
        co_parent_coverage.append(co_parent_cov[pos])
        progeny_coverage.append(progeny_cov[pos])

    timers.stop("arrangement")

    # create plots
    if param.plot:
        timers.start("plots")
        create_genetic_plots(
            param.output,
            chromosome,
            sample,
            ground_truth_reader,
            varinfo,
            clustering,
            haplo_skeletons,
            haplotypes,
            phased_positions,
            parent_coverage,
            co_parent_coverage,
            progeny_coverage,
            param,
        )
        timers.stop("plots")

    return superreads, components


def determine_pedigree(pedigree_file, samples, parent_samples, progeny_samples=None):
    """
    Evaluates the pedigree file and returns:
    1.  The final set of samples to phase (Only different from input when empty)
    2.  Dictionary that maps samples to co-parents
    3.  Dictionary that maps samples to list of progenies

    Progenies are only accounted for when defined as trio in pedigree file and present in
    progeny VCF if present or in arental/primary VCF otherwise.

    pedigree_file -- Input pedigree file
    samples -- User defined samples to phase. If None, will be inferred from input
    parent_samples -- List of samples in parental/primary VCF
    progeny_samples -- List of samples in progeny VCF. None if no progeny VCF given
    """

    # For each parent, store co-parent and list of progenies
    coparents = dict()
    progenies = dict()
    with open(pedigree_file, "r") as ped:
        for i, line in enumerate(ped):
            tokens = line.replace("\n", "").split(" ")
            if len(tokens) != 3:
                logger.error(f"Line {i} in pedfile contains {len(tokens)} values instead of 3.")
                raise CommandLineError(None)
            # parent = tokens[0]
            # co_parent = tokens[1]
            progeny = tokens[2]
            if progeny in tokens[:2]:
                logger.warning(f"Ignore: Sample {progeny} defined as its own parent in line {i}.")
                continue
            for parent, co_parent in zip(tokens[:2], tokens[-2::-1]):
                if parent in coparents and coparents[parent] != co_parent:
                    other = coparents[parent]
                    msg = f"Pedfile assigns multiple partners ({co_parent}, {other}) to {parent}. Currently only one partner per sample is supported."
                    logger.error(msg)
                    raise CommandLineError(msg)
                coparents[parent] = co_parent
                if parent not in progenies:
                    progenies[parent] = []
                if progeny in progenies[parent]:
                    logger.warning(
                        f"Ignore: Duplicate trio ({parent}, {co_parent}, {progeny}) in pedfile line {i}"
                    )
                else:
                    progenies[parent].append(progeny)

    # Validate:
    # 1: Each requested phasable sample must occur as parent in pediegree file
    # 2: Each parent must have exactly one co-parent occuring in pedigree file AND parental VCF
    if samples:
        # Check every requested sample to have co-parent and both be defined in parent VCF
        for sample in samples:
            if sample not in coparents:
                msg = f"Requested parent sample {sample} does not occur as parent in pedfile."
                logger.error(msg)
                raise CommandLineError(msg)
            if sample not in parent_samples:
                msg = f"Requested parent sample {sample} is not present in primary VCF file."
                logger.error(msg)
                raise CommandLineError(msg)
            if coparents[parent] not in parent_samples:
                msg = f"Partner {coparents[parent]} of requested parent sample {sample} is not present in primary VCF file."
                logger.error(msg)
                raise CommandLineError(msg)
    else:
        # If not specified, find all phasable parents. Raise error if no sample applicable
        samples = []
        if not coparents:
            msg = "Pedfile does not contain any trios."
            logger.error(msg)
            raise CommandLineError(msg)
        for sample in coparents:
            if sample in parent_samples:
                samples.append(sample)
        if not samples:
            msg = "No prospect parent sample from the pedfile is present in primary VCF file"
            logger.error(msg)
            raise CommandLineError(msg)

    # filter out progeny samples unpresent in VCFs and unrequested parent samples
    fprogenies = dict()
    fcoparents = dict()
    for sample in samples:
        fprogenies[sample] = []
        fcoparents[sample] = coparents[sample]
        for progeny in progenies[sample]:
            # if progeny file given: use only this and ignore progeny samples from primary VCF
            if progeny_samples:
                if progeny in progeny_samples:
                    fprogenies[sample].append(progeny)
                elif progeny in parent_samples:
                    logger.warning(
                        f"Ignore: Progeny {progeny} present in primary VCF instead of progeny VCF."
                    )
                else:
                    logger.warning(f"Ignore: Progeny {progeny} not present in progeny VCF.")
            else:
                if progeny in parent_samples:
                    fprogenies[sample].append(progeny)
                else:
                    logger.warning(f"Ignore: Progeny {progeny} not present in primary VCF.")

    return samples, fcoparents, fprogenies


def get_parent_progeny_coverage(parent, co_parent, progeny_list, parent_table, progeny_table):
    # store progeny coverage
    parent_depths = parent_table.allele_depths_of(parent)
    co_parent_depths = parent_table.allele_depths_of(co_parent)
    parent_cov = [sum(parent_depths[pos]) for pos in range(len(parent_table))]
    co_parent_cov = [sum(co_parent_depths[pos]) for pos in range(len(parent_table))]
    progeny_cov = [0 for _ in range(len(parent_table))]
    for off in progeny_list:
        parent_pos = 0
        progeny_pos = 0
        allele_depths = progeny_table.allele_depths_of(off)
        assert len(allele_depths) == len(progeny_table)
        while progeny_pos < len(allele_depths) and parent_pos < len(parent_table):
            if (
                parent_table.variants[parent_pos].position
                == progeny_table.variants[progeny_pos].position
            ):
                progeny_cov[parent_pos] += sum(allele_depths[progeny_pos])
                progeny_pos += 1
            else:
                assert (
                    parent_table.variants[parent_pos].position
                    < progeny_table.variants[progeny_pos].position
                )
            parent_pos += 1
    return parent_cov, co_parent_cov, progeny_cov


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
    arg("--indels", dest="indels_used", action="store_true", help=argparse.SUPPRESS)
    arg("--only-snvs", action="store_true", help="Phase only SNVs")
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
        "--scoring-window",
        metavar="SCORINGWINDOW",
        dest="scoring_window",
        type=int,
        default=250,
        required=False,
        help="Size of the window (in variants) for statistical progeny scoring.",
    )
    arg(
        "--complexity-support",
        "-C",
        dest="complexity_support",
        type=int,
        default=0,
        required=False,
        help="Indicates what level of genotype complexity is allowed for phased variants. 0 = simplex-nulliplex only, 1 = simplex-simplex on top, 2 = duplex-nulliplex on top. Default is 0.",
    )
    arg(
        "--distrust-genotypes",
        dest="distrust_genotypes",
        default=False,
        action="store_true",
        help="Internally retypes the reported parent genotypes based on allele distribution in progeny samples.",
    )

    # more arguments, which are experimental or for debugging and should not be presented to the user
    arg(
        "--ratio-cutoff",
        metavar="RATIOCUTOFF",
        dest="ratio_cutoff",
        type=float,
        default=0.0,
        required=False,
        help=argparse.SUPPRESS,
    )  # help="Cutoff threshold for deviation of coverage ratio between parent, co-parent and progeny compared to median."
    arg(
        "--allele-error-rate",
        metavar="ALLELEERRORRATE",
        dest="allele_error_rate",
        type=float,
        default=0.06,
        required=False,
        help=argparse.SUPPRESS,
    )  # help="Assumed error rate for observed alleles among parents and progeny."
    arg(
        "--plot",
        dest="plot",
        default=False,
        action="store_true",
        help=argparse.SUPPRESS,
    )  # help="Plots the variant clustering and arrangement."
    arg(
        "--ground-truth-file",
        "-g",
        required=False,
        help=argparse.SUPPRESS,
    )  # help="File containing a ground truth phasing. Only used for plotting."


def validate(args, parser):
    if args.allele_error_rate > 0.5 or args.allele_error_rate < 0.01:
        parser.error("Allele error rate must be between 0.01 and 0.5.")
    if args.scoring_window < 1:
        parser.error("Scoring window must be a positive integer.")
    if args.complexity_support not in [0, 1, 2]:
        parser.error("Complexity support level must be either 0, 1 or 2.")
    if args.ploidy % 2 > 0:
        parser.error("Odd ploidies are not supported.")
    if args.ploidy < 2:
        parser.error("Ploidy must be at least 2.")
    if args.indels_used:
        logger.warning("Ignoring --indels as indel phasing is default in WhatsHap 2.0+")


def main(args):
    run_polyphasegenetic(**vars(args))
