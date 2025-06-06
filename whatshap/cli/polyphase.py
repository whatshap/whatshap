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

from contextlib import ExitStack

from typing import (
    Optional,
    List,
    TextIO,
    Sequence,
    FrozenSet,
    Dict,
    Tuple,
)

from whatshap import __version__
from whatshap.core import (
    Read,
    ReadSet,
    NumericSampleIds,
)
from whatshap.cli import log_memory_usage, PhasedInputReader, CommandLineError

from whatshap.polyphase import (
    PolyphaseParameter,
    create_genotype_list,
    extract_partial_phasing,
    Position,
)
from whatshap.polyphase.algorithm import solve_polyphase_instance, compute_cut_positions
from whatshap.polyphase.plots import draw_plots
from whatshap.polyphase.solver import AlleleMatrix

from whatshap.timer import StageTimer
from whatshap.utils import ChromosomeFilter
from whatshap.vcf import VcfReader, PhasedVcfWriter, VariantTable, PloidyError

__author__ = "Jana Ebler, Sven Schrinner"


logger = logging.getLogger(__name__)


def run_polyphase(
    phase_input_files: Sequence[str],
    variant_file: str,
    ploidy: int,
    reference: Optional[str] = None,
    output: TextIO = sys.stdout,
    samples: Optional[Sequence[str]] = None,
    chromosomes: Optional[List[str]] = None,
    excluded_chromosomes: Optional[List[str]] = None,
    ignore_read_groups: bool = False,
    only_snvs: bool = False,
    mapping_quality: int = 20,
    distrust_genotypes: bool = False,
    tag: str = "PS",
    read_list_filename: Optional[str] = None,
    write_command_line_header: bool = True,
    include_haploid_sets: bool = False,
    block_cut_sensitivity=4,
    min_overlap: int = 2,
    mav: bool = True,
    threads: int = 1,
    use_supplementary: bool = False,
    supplementary_distance_threshold: int = 100_000,
    use_prephasing: bool = False,
    ce_bundle_edges: bool = False,
    plot_clusters: bool = False,
    plot_threading: bool = False,
):
    """
    Run Polyploid Phasing.

    phase_input_files -- list of paths to BAM/CRAM/VCF files
    variant-file -- path to input VCF
    ploidy -- target ploidy for all chromosomes
    reference -- path to reference FASTA
    output -- path to output VCF or a file like object
    samples -- names of samples to phase. An empty list means: phase all samples
    chromosomes -- names of chromosomes to phase. An empty list means: phase all chromosomes
    excluded_chromosomes -- names of chromosomes not to phase.
    ignore_read_groups -- assigns all reads to one sample. Cannot be used for multi-sample VCF
    mapping_quality -- discard reads below this mapping quality
    distrust_genotypes -- allows to override the provided genotypes from VCF in phasing output
    tag -- How to store phasing info in the VCF, can be 'PS' or 'HP'
    read_list_filename -- name of file to write list of used reads to
    write_command_line_header -- whether to add a ##commandline header to the output VCF
    include_haploid_sets -- writes independent phase block identifiers for each phase
    block_cut_sensitivity -- policy for creating block cuts
    min_overlap -- minimum number of common variants to of read pair for score computation
    mav -- include multi-allelic variants
    threads -- number of worker threads to process disconnected blocks
    use_supplementary -- use supplementary alignments with primary
    supplementary_distance_threshold -- distance threshold for filtering supplementary alignments
    use_prephasing -- consider existing phasing in input VCF
    ce_bundle_edges -- alternative edge contraction policy in cluster editing heuristic
    plot_clusters -- add plot of cluster editing result
    plot_threading -- add plot for threading result
    verify_genotypes -- posterior genotype correction based on observed allele frequency
    """
    timers = StageTimer()
    logger.info(
        "This is WhatsHap (polyploid) %s running under Python %s",
        __version__,
        platform.python_version(),
    )
    numeric_sample_ids = NumericSampleIds()
    with ExitStack() as stack:
        assert phase_input_files
        phased_input_reader = stack.enter_context(
            PhasedInputReader(
                phase_input_files,
                reference,
                numeric_sample_ids,
                ignore_read_groups,
                only_snvs=only_snvs,
                mapq_threshold=mapping_quality,
                use_supplementary=use_supplementary,
                supplementary_distance_threshold=supplementary_distance_threshold,
            )
        )
        assert not phased_input_reader.has_vcfs

        if write_command_line_header:
            command_line = "(whatshap {}) {}".format(__version__, " ".join(sys.argv[1:]))
        else:
            command_line = None
        try:
            vcf_writer: PhasedVcfWriter = stack.enter_context(
                PhasedVcfWriter(
                    command_line=command_line,
                    in_path=variant_file,
                    out_file=output,
                    tag=tag,
                    ploidy=ploidy,
                    only_snvs=only_snvs,
                    include_haploid_sets=include_haploid_sets,
                    mav=mav,
                )
            )
        except OSError as e:
            raise CommandLineError(e)

        vcf_reader = stack.enter_context(
            VcfReader(
                variant_file,
                only_snvs=only_snvs,
                phases=True,
                genotype_likelihoods=False,
                ploidy=ploidy,
                mav=mav,
            )
        )

        if ignore_read_groups and not samples and len(vcf_reader.samples) > 1:
            raise CommandLineError(
                "When using --ignore-read-groups on a VCF with "
                "multiple samples, --sample must also be used."
            )
        if not samples:
            samples = vcf_reader.samples

        vcf_sample_set = set(vcf_reader.samples)
        for sample in samples:
            if sample not in vcf_sample_set:
                raise CommandLineError(
                    f"Sample {sample!r} requested on command-line not found in VCF"
                )

        if use_prephasing and block_cut_sensitivity > 1:
            logger.info(
                "Consider using '-B 0' or '-B 1' when adding pre-phasings from another source."
            )

        samples = frozenset(samples)

        read_list_file = None
        if read_list_filename:
            raise NotImplementedError("create_read_list_file not implemented")
            # read_list_file = create_read_list_file(read_list_filename)

        # Store phasing parameters in tuple to keep function signatures cleaner
        phasing_param = PolyphaseParameter(
            ploidy=ploidy,
            ce_bundle_edges=ce_bundle_edges,
            distrust_genotypes=distrust_genotypes,
            min_overlap=min_overlap,
            block_cut_sensitivity=block_cut_sensitivity,
            plot_clusters=plot_clusters,
            plot_threading=plot_threading,
            plot_path=output if type(output) is str else output.name,
            threads=threads,
            use_prephasing=use_prephasing,
        )

        try:
            included_chromosomes = ChromosomeFilter(chromosomes, excluded_chromosomes)
            for variant_table in timers.iterate("parse_vcf", vcf_reader):
                chromosome = variant_table.chromosome
                if chromosome in included_chromosomes:
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

                # These three variables hold the phasing results for all samples
                components, haploid_components, superreads = phase_single_chromosome(
                    variant_table, phased_input_reader, samples, timers, phasing_param
                )
                # Unphasable variants are removed from input table from here!

                with timers("write_vcf"):
                    logger.info("======== Writing VCF")
                    vcf_writer.write(
                        chromosome,
                        superreads,
                        components,
                        haploid_components if include_haploid_sets else None,
                    )
                    logger.info("Done writing VCF")
                logger.debug("Chromosome %r finished", chromosome)
        except PloidyError as e:
            raise CommandLineError(e)

    if read_list_file:
        read_list_file.close()

    logger.info("\n== SUMMARY ==")

    log_memory_usage(include_children=(threads > 1))
    logger.info("Time spent reading BAM/CRAM:         %6.1f s", timers.elapsed("read_bam"))
    logger.info("Time spent parsing VCF:              %6.1f s", timers.elapsed("parse_vcf"))
    logger.info("Time spent detecting blocks:         %6.1f s", timers.elapsed("detecting_blocks"))
    if threads == 1:
        logger.info("Time spent scoring reads:            %6.1f s", timers.elapsed("read_scoring"))
        logger.info("Time spent solving cluster editing:  %6.1f s", timers.elapsed("clustering"))
        logger.info("Time spent threading haplotypes:     %6.1f s", timers.elapsed("threading"))
        logger.info("Time spent reordering haplotypes:    %6.1f s", timers.elapsed("reordering"))
    else:
        logger.info("Time spent phasing blocks:           %6.1f s", timers.elapsed("phase_blocks"))
    if plot_clusters or plot_threading:
        logger.info("Time spent creating plots:           %6.1f s", timers.elapsed("create_plots"))
    logger.info("Time spent writing VCF:              %6.1f s", timers.elapsed("write_vcf"))
    logger.info("Time spent on rest:                  %6.1f s", timers.total() - timers.sum())
    logger.info("Total elapsed time:                  %6.1f s", timers.total())


def phase_single_chromosome(
    variant_table: VariantTable,
    phased_input_reader: PhasedInputReader,
    samples: FrozenSet[str],
    timers: StageTimer,
    param: PolyphaseParameter,
) -> Tuple[
    Dict[str, Dict[Position, int]], Dict[str, Dict[Position, List[int]]], Dict[str, ReadSet]
]:
    chromosome = variant_table.chromosome
    superreads, components, haploid_components = dict(), dict(), dict()

    # Iterate over all samples to process
    for sample in samples:
        logger.info("---- Processing individual %s", sample)

        # Process inputs for this sample
        missing_genotypes = set()
        heterozygous = set()

        genotypes = variant_table.genotypes_of(sample)
        for index, gt in enumerate(genotypes):
            if gt.is_none():
                missing_genotypes.add(index)
            elif not gt.is_homozygous():
                heterozygous.add(index)
            else:
                assert gt.is_homozygous()
        to_discard = set(range(len(variant_table))).difference(heterozygous)
        # Remove calls to be discarded from variant table
        phasable_variant_table = variant_table.create_subtable([sample])
        phasable_variant_table.remove_rows_by_index(to_discard)

        logger.info(
            "Number of variants skipped due to missing genotypes: %d",
            len(missing_genotypes),
        )
        logger.info("Number of remaining heterozygous variants: %d", len(phasable_variant_table))

        if len(phasable_variant_table) < 2:
            logger.debug("Skipped phasing because there is only one variant")
            continue

        # Get the reads belonging to this sample
        timers.start("read_bam")
        readset, vcf_source_ids = phased_input_reader.read(
            chromosome, phasable_variant_table.variants, sample
        )
        readset.sort()
        timers.stop("read_bam")

        # Remove reads with insufficient variants
        readset = readset.subset(
            [i for i, read in enumerate(readset) if len(read) >= max(2, param.min_overlap)]
        )
        if len(readset) == 0:
            logger.debug("Skipped phasing because no suitable reads remain")
            continue

        logger.info("Kept %d reads that cover at least two variants each", len(readset))

        # Adapt the variant table to the subset of reads
        phasable_variant_table.subset_rows_by_position(readset.get_positions())

        # Run the actual phasing
        (
            sample_components,
            sample_haploid_components,
            sample_superreads,
        ) = phase_single_individual(readset, phasable_variant_table, sample, param, timers)
        del phasable_variant_table

        # Collect results
        components[sample] = sample_components
        haploid_components[sample] = sample_haploid_components
        superreads[sample] = sample_superreads

    return components, haploid_components, superreads


def phase_single_individual(
    readset: ReadSet,
    phasable_variant_table: VariantTable,
    sample: str,
    param: PolyphaseParameter,
    timers: StageTimer,
) -> Tuple[Dict[Position, int], Dict[Position, List[int]], ReadSet]:
    # Compute the genotypes that belong to the variant table and create a list of all genotypes
    genotype_list = create_genotype_list(phasable_variant_table, sample)

    # Optional: Extract partial phasing from variant table
    prephasing = None
    if param.use_prephasing:
        prephasing = extract_partial_phasing(phasable_variant_table, sample, param.ploidy)
        if prephasing is None:
            logger.warning(
                f"Input VCF does not contain any phased blocks for {sample}. "
                "No pre-phasing will be used for this sample."
            )

    # Retrieve solution
    allele_matrix = AlleleMatrix(readset)
    result = solve_polyphase_instance(allele_matrix, genotype_list, param, timers, prephasing)
    cuts, hap_cuts = compute_cut_positions(
        result.breakpoints, param.ploidy, param.block_cut_sensitivity
    )

    # Summarize data for VCF file
    accessible_pos = sorted(readset.get_positions())
    components = {}
    haploid_components = {}

    num_vars = len(readset.get_positions())
    cuts = cuts + [num_vars]
    for i, cut_pos in enumerate(cuts[:-1]):
        for pos in range(cuts[i], cuts[i + 1]):
            components[accessible_pos[pos]] = accessible_pos[cuts[i]]
            components[accessible_pos[pos] + 1] = accessible_pos[cuts[i]]
            haploid_components[accessible_pos[pos]] = [0] * param.ploidy
            haploid_components[accessible_pos[pos] + 1] = [0] * param.ploidy

    for j in range(param.ploidy):
        hap_cuts[j] = hap_cuts[j] + [num_vars]
        for i, cut_pos in enumerate(hap_cuts[j][:-1]):
            for pos in range(hap_cuts[j][i], hap_cuts[j][i + 1]):
                haploid_components[accessible_pos[pos]][j] = accessible_pos[hap_cuts[j][i]]
                haploid_components[accessible_pos[pos] + 1][j] = accessible_pos[hap_cuts[j][i]]

    # insert alleles
    superreads = ReadSet()
    phased_pos = [i for i in range(num_vars) if -1 not in [h[i] for h in result.haplotypes]]
    for i in range(param.ploidy):
        read = Read(f"superread {i + 1}", 0, 0)
        for j in phased_pos:
            read.add_variant(accessible_pos[j], result.haplotypes[i][j], 0)
        superreads.add(read)

    # Plot option
    if param.plot_clusters or param.plot_threading:
        timers.start("create_plots")
        draw_plots(
            readset,
            result,
            cuts[:-1],
            phasable_variant_table,
            param.plot_clusters,
            param.plot_threading,
            param.plot_path,
        )
        timers.stop("create_plots")

    # Return results
    return components, haploid_components, superreads


def add_arguments(parser):
    arg = parser.add_argument
    # Positional argument
    arg(
        "variant_file",
        metavar="VCF",
        help="VCF file with variants to be phased (can be gzip-compressed)",
    )
    arg(
        "phase_input_files",
        nargs="*",
        metavar="PHASEINPUT",
        help="BAM or CRAM with sequencing reads.",
    )
    arg(
        "-o",
        "--output",
        default=sys.stdout,
        help="Output VCF file. Add .gz to the file name to get compressed output. "
        "If omitted, use standard output.",
    )
    arg(
        "--reference",
        "-r",
        metavar="FASTA",
        help="Reference file. Provide this to detect alleles through re-alignment. "
        "If no index (.fai) exists, it will be created",
    )
    arg(
        "--tag",
        choices=("PS", "HP"),
        default="PS",
        help="Store phasing information with PS tag (standardized) or "
        "HP tag (used by GATK ReadBackedPhasing) (default: %(default)s)",
    )
    # arg(
    #    "--output-read-list",
    #    metavar="FILE",
    #    default=None,
    #    dest="read_list_filename",
    #    help="Write reads that have been used for phasing to FILE.",
    # )

    arg = parser.add_argument_group("Input pre-processing, selection, and filtering").add_argument
    arg(
        "--mapping-quality",
        "--mapq",
        metavar="QUAL",
        default=20,
        type=int,
        help="Minimum mapping quality (default: %(default)s)",
    )
    arg("--indels", dest="indels_used", action="store_true", help=argparse.SUPPRESS)
    arg("--only-snvs", action="store_true", help="Only phase SNVs")
    arg(
        "--ignore-read-groups",
        default=False,
        action="store_true",
        help="Ignore read groups in BAM/CRAM header and assume all reads come "
        "from the same sample.",
    )
    arg(
        "--include-haploid-sets",
        default=False,
        action="store_true",
        help="Include the phase set information for every single haplotype in a custom VCF format field 'HS'.",
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
    arg(
        "--exclude-chromosome",
        dest="excluded_chromosomes",
        default=[],
        action="append",
        help="Name of chromosome not to phase.",
    )
    arg(
        "--distrust-genotypes",
        dest="distrust_genotypes",
        action="store_true",
        default=False,
        help="Allows the phaser to change genotypes if beneficial for the internal model.",
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
        "--use-prephasing",
        dest="use_prephasing",
        action="store_true",
        default=False,
        help="Uses existing phase set blocks in the input to increase contiguity of phasing output.",
    )
    arg(
        "--min-overlap",
        metavar="OVERLAP",
        type=int,
        default=2,
        help="Minimum required read overlap for internal read clustering stage (default: %(default)s).",
    )
    arg(
        "--block-cut-sensitivity",
        "-B",
        metavar="SENSITIVITY",
        type=int,
        dest="block_cut_sensitivity",
        default=4,
        help="Strategy to determine block borders. 0 yields the longest blocks with more switch errors, 5 has the shortest blocks with lowest switch error rate (default: %(default)s).",
    )
    arg(
        "--threads",
        "-t",
        metavar="THREADS",
        type=int,
        default=1,
        help="Maximum number of CPU threads used (default: %(default)s).",
    )
    arg(
        "--no-mav",
        dest="mav",
        default=True,
        action="store_false",
        help="Disables phasing of multi-allelic variants.",
    )
    arg(
        "--use-supplementary",
        dest="use_supplementary",
        action="store_true",
        default=False,
        help="Use also supplementary alignments (default: ignore supplementary_ alignments)",
    )
    arg(
        "--supplementary-distance",
        metavar="DIST",
        type=int,
        dest="supplementary_distance_threshold",
        default=100_000,
        help="Skip supplementary alignments further than DIST bp away from the primary alignment (default: %(default)s)",
    )

    # more arguments, which are experimental or for debugging and should not be presented to the user
    arg(
        "--ce-bundle-edges",
        dest="ce_bundle_edges",
        default=False,
        action="store_true",
        help=argparse.SUPPRESS,
    )  # help='Influences the cluster editing heuristic. Only for debug/developing purpose (default: %(default)s).')
    arg(
        "--plot-clusters",
        dest="plot_clusters",
        default=False,
        action="store_true",
        help=argparse.SUPPRESS,
    )  # help='Plot a super heatmap for the computed clustering (default: %(default)s).')
    arg(
        "--plot-threading",
        dest="plot_threading",
        default=False,
        action="store_true",
        help=argparse.SUPPRESS,
    )  # help='Plot the haplotypes\' threading through the read clusters (default: %(default)s).')
    arg(
        "--verify-genotypes",
        default=False,
        action="store_true",
        help=argparse.SUPPRESS,
    )  # help="Verify input genotypes by re-typing them using the given reads.",


def validate(args, parser):
    if args.block_cut_sensitivity > 5 or args.block_cut_sensitivity < 0:
        parser.error("Block cut sensitivity must be an integer value between 0 and 5.")
    if args.indels_used:
        logger.warning("Ignoring --indels as indel phasing is default in WhatsHap 2.0+")
    if args.verify_genotypes:
        logger.warning("Ignoring deprecated option --verify-genotypes.")


def main(args):
    del args.indels_used
    del args.verify_genotypes
    run_polyphase(**vars(args))
