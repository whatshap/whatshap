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

from collections import namedtuple
from copy import deepcopy
from multiprocessing import Pool

from contextlib import ExitStack

from whatshap import __version__
from whatshap.core import (
    Read,
    ReadSet,
    Genotype,
    ClusterEditingSolver,
    NumericSampleIds,
    compute_polyploid_genotypes,
    scoreReadset,
)
from whatshap.cli import log_memory_usage, PhasedInputReader, CommandLineError
from whatshap.polyphaseplots import draw_plots
from whatshap.polyphaseutil import create_genotype_list, compute_block_starts, split_readset
from whatshap.reorder import run_reordering
from whatshap.threading import run_threading
from whatshap.timer import StageTimer
from whatshap.vcf import VcfReader, PhasedVcfWriter, PloidyError

__author__ = "Jana Ebler, Sven Schrinner"

PhasingParameter = namedtuple(
    "PhasingParameter",
    [
        "ploidy",
        "verify_genotypes",
        "ce_bundle_edges",
        "distrust_genotypes",
        "min_overlap",
        "block_cut_sensitivity",
        "plot_clusters",
        "plot_threading",
        "threads",
    ],
)

SingleBlockResult = namedtuple(
    "SingleBlockResult",
    ["block_id", "clustering", "paths", "cuts", "hap_cuts", "haplotypes"],
)

logger = logging.getLogger(__name__)


def run_polyphase(
    phase_input_files,
    variant_file,
    ploidy,
    reference=None,
    output=sys.stdout,
    samples=None,
    chromosomes=None,
    verify_genotypes=False,
    ignore_read_groups=False,
    indels=True,
    mapping_quality=20,
    tag="PS",
    include_haploid_sets=False,
    distrust_genotypes=False,
    write_command_line_header=True,
    read_list_filename=None,
    ce_bundle_edges=False,
    min_overlap=2,
    plot_clusters=False,
    plot_threading=False,
    block_cut_sensitivity=4,
    threads=1,
):
    """
    Run Polyploid Phasing.

    phase_input_files -- list of paths to BAM/CRAM/VCF files
    variant-file -- path to input VCF
    reference -- path to reference FASTA
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
                indels=indels,
                mapq_threshold=mapping_quality,
            )
        )
        assert not phased_input_reader.has_vcfs

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
                    include_haploid_sets=include_haploid_sets,
                )
            )
        except OSError as e:
            raise CommandLineError(e)

        vcf_reader = stack.enter_context(
            VcfReader(
                variant_file, indels=indels, phases=True, genotype_likelihoods=False, ploidy=ploidy
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
                    "Sample {!r} requested on command-line not found in VCF".format(sample)
                )

        samples = frozenset(samples)

        read_list_file = None
        if read_list_filename:
            raise NotImplementedError("create_read_list_file not implemented")
            # read_list_file = create_read_list_file(read_list_filename)

        # Store phasing parameters in tuple to keep function signatures cleaner
        phasing_param = PhasingParameter(
            ploidy=ploidy,
            verify_genotypes=verify_genotypes,
            ce_bundle_edges=ce_bundle_edges,
            distrust_genotypes=distrust_genotypes,
            min_overlap=min_overlap,
            block_cut_sensitivity=block_cut_sensitivity,
            plot_clusters=plot_clusters,
            plot_threading=plot_threading,
            threads=threads,
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
                    phasable_variant_table = deepcopy(variant_table)
                    # Remove calls to be discarded from variant table
                    phasable_variant_table.remove_rows_by_index(to_discard)

                    logger.info(
                        "Number of variants skipped due to missing genotypes: %d",
                        len(missing_genotypes),
                    )
                    logger.info(
                        "Number of remaining heterozygous variants: %d", len(phasable_variant_table)
                    )

                    # Get the reads belonging to this sample
                    timers.start("read_bam")
                    readset, vcf_source_ids = phased_input_reader.read(
                        chromosome, phasable_variant_table.variants, sample
                    )
                    readset.sort()
                    timers.stop("read_bam")

                    # Verify genotypes
                    if verify_genotypes:
                        timers.start("verify_genotypes")
                        logger.info("Verify genotyping of %s", sample)
                        positions = [v.position for v in phasable_variant_table.variants]
                        computed_genotypes = [
                            Genotype(gt)
                            for gt in compute_polyploid_genotypes(readset, ploidy, positions)
                        ]
                        # skip all positions at which genotypes do not match
                        given_genotypes = phasable_variant_table.genotypes_of(sample)
                        matching_genotypes = []
                        missing_genotypes = set()
                        for i, g in enumerate(given_genotypes):
                            c_g = computed_genotypes[i]
                            if (g == c_g) or (c_g is None):
                                matching_genotypes.append(g)
                            else:
                                matching_genotypes.append(Genotype([]))
                                missing_genotypes.add(i)
                        phasable_variant_table.set_genotypes_of(sample, matching_genotypes)

                        # Remove variants with deleted genotype
                        phasable_variant_table.remove_rows_by_index(missing_genotypes)
                        logger.info(
                            "Number of variants removed due to inconsistent genotypes: %d",
                            len(missing_genotypes),
                        )
                        logger.info(
                            "Number of remaining heterozygous variants: %d",
                            len(phasable_variant_table),
                        )

                        # Re-read the readset to remove discarded variants
                        readset, vcf_source_ids = phased_input_reader.read(
                            chromosome, phasable_variant_table.variants, sample
                        )
                        readset.sort()
                        timers.stop("verify_genotypes")

                    # Remove reads with insufficient variants
                    readset = readset.subset(
                        [i for i, read in enumerate(readset) if len(read) >= max(2, min_overlap)]
                    )
                    logger.info("Kept %d reads that cover at least two variants each", len(readset))

                    # Adapt the variant table to the subset of reads
                    phasable_variant_table.subset_rows_by_position(readset.get_positions())

                    # Run the actual phasing
                    (
                        sample_components,
                        sample_haploid_components,
                        sample_superreads,
                    ) = phase_single_individual(
                        readset, phasable_variant_table, sample, phasing_param, output, timers
                    )

                    # Collect results
                    components[sample] = sample_components
                    haploid_components[sample] = sample_haploid_components
                    superreads[sample] = sample_superreads

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
                timers.start("parse_vcf")
            timers.stop("parse_vcf")
        except PloidyError as e:
            raise CommandLineError(e)

    if read_list_file:
        read_list_file.close()

    logger.info("\n== SUMMARY ==")

    log_memory_usage(include_children=(threads > 1))
    logger.info("Time spent reading BAM/CRAM:         %6.1f s", timers.elapsed("read_bam"))
    logger.info("Time spent parsing VCF:              %6.1f s", timers.elapsed("parse_vcf"))
    if verify_genotypes:
        logger.info(
            "Time spent verifying genotypes:      %6.1f s", timers.elapsed("verify_genotypes")
        )
    logger.info("Time spent detecting blocks:         %6.1f s", timers.elapsed("detecting_blocks"))
    if threads == 1:
        logger.info("Time spent scoring reads:            %6.1f s", timers.elapsed("read_scoring"))
        logger.info("Time spent solving cluster editing:  %6.1f s", timers.elapsed("clustering"))
        logger.info("Time spent threading haplotypes:     %6.1f s", timers.elapsed("threading"))
        logger.info("Time spent reordering haplotypes:    %6.1f s", timers.elapsed("reordering"))
    else:
        """
        TODO: The runtime measurement for the different stages does not properly for multithreading,
        because the global timer is not visible from within the phase_single_block_mt method.
        Workaround is to only report the total phasing time.
        """
        logger.info("Time spent phasing blocks:           %6.1f s", timers.elapsed("phase_blocks"))
    if plot_clusters or plot_threading:
        logger.info("Time spent creating plots:           %6.1f s", timers.elapsed("create_plots"))
    logger.info("Time spent writing VCF:              %6.1f s", timers.elapsed("write_vcf"))
    logger.info("Time spent on rest:                  %6.1f s", timers.total() - timers.sum())
    logger.info("Total elapsed time:                  %6.1f s", timers.total())


def phase_single_individual(readset, phasable_variant_table, sample, param, output, timers):

    # Compute the genotypes that belong to the variant table and create a list of all genotypes
    genotype_list = create_genotype_list(phasable_variant_table, sample)

    # Precompute block borders based on read coverage and linkage between variants
    logger.info("Detecting connected components with weak interconnect ..")
    timers.start("detecting_blocks")
    num_vars = len(readset.get_positions())
    ploidy = param.ploidy
    if param.block_cut_sensitivity == 0:
        block_starts = [0]
    elif param.block_cut_sensitivity == 1:
        block_starts = compute_block_starts(readset, ploidy, single_linkage=True)
    else:
        block_starts = compute_block_starts(readset, ploidy, single_linkage=False)

    # Set block borders and split readset
    block_starts.append(num_vars)
    num_non_singleton_blocks = len(
        [i for i in range(len(block_starts) - 1) if block_starts[i] < block_starts[i + 1] - 1]
    )
    logger.info(
        "Split heterozygous variants into {} blocks (and {} singleton blocks).".format(
            num_non_singleton_blocks, len(block_starts) - num_non_singleton_blocks - 1
        )
    )

    block_readsets = split_readset(readset, block_starts)
    timers.stop("detecting_blocks")

    # Process blocks independently
    results = []

    # Create genotype slices for blocks
    gt_slices = []
    for block_id, block_readset in enumerate(block_readsets):
        block_start = block_starts[block_id]
        block_end = block_starts[block_id + 1]
        block_num_vars = block_end - block_start

        assert len(block_readset.get_positions()) == block_num_vars
        gt_slices.append(genotype_list[block_start:block_end])

    processed_non_singleton_blocks = 0

    """
    Python's multiprocessing makes hard copies of the passed arguments, which is not trivial for
    cython objects, especially when they contain pointers to other cython objects. Any passed
    object must be (de)serializable (in Python: pickle). All other objects created in the main
    thread are also accessible by the workers, but they are handled via the copy-on-write policy.
    This means, that e.g. the large main readset is not hardcopied for every thread, as long as it
    is not modified there. Since this would cause a massive waste of memory, this must not be done
    and the main readset must also never be passed as argument to the workers.
    """
    if param.threads == 1:
        # for single-threading, process everything individually to minimize memory footprint
        for block_id, block_readset in enumerate(block_readsets):
            block_num_vars = block_starts[block_id + 1] - block_starts[block_id]
            if block_num_vars > 1:
                processed_non_singleton_blocks += 1
                logger.info(
                    "Processing block {} of {} with {} reads and {} variants.".format(
                        processed_non_singleton_blocks,
                        num_non_singleton_blocks,
                        len(block_readset),
                        block_num_vars,
                    )
                )
            results.append(
                phase_single_block(block_id, block_readset, gt_slices[block_id], param, timers)
            )

    else:
        # sort block by descending size (4/3-approximation for scheduling problem)
        joblist = [(i, len(block_readsets[i])) for i in range(len(block_readsets))]
        joblist.sort(key=lambda x: -x[1])

        timers.start("phase_blocks")
        with Pool(processes=param.threads) as pool:
            process_results = [
                pool.apply_async(
                    phase_single_block_mt,
                    (
                        block_id,
                        block_readsets[block_id],
                        gt_slices[block_id],
                        param,
                        timers,
                        job_id,
                        num_non_singleton_blocks,
                    ),
                )
                for job_id, (block_id, block_readset) in enumerate(joblist)
            ]
            # collect all blockwise results
            blockwise_results = [res.get() for res in process_results]
            results = sorted(blockwise_results, key=lambda x: x.block_id)

        timers.stop("phase_blocks")

    # Aggregate blockwise results
    clustering, threading, haplotypes, cuts, hap_cuts = aggregate_results(results, ploidy)

    # Summarize data for VCF file
    accessible_pos = sorted(readset.get_positions())
    components = {}
    haploid_components = {}

    cuts = cuts + [num_vars]
    for i, cut_pos in enumerate(cuts[:-1]):
        for pos in range(cuts[i], cuts[i + 1]):
            components[accessible_pos[pos]] = accessible_pos[cuts[i]]
            components[accessible_pos[pos] + 1] = accessible_pos[cuts[i]]
            haploid_components[accessible_pos[pos]] = [0] * ploidy
            haploid_components[accessible_pos[pos] + 1] = [0] * ploidy

    for j in range(ploidy):
        hap_cuts[j] = hap_cuts[j] + [num_vars]
        for i, cut_pos in enumerate(hap_cuts[j][:-1]):
            for pos in range(hap_cuts[j][i], hap_cuts[j][i + 1]):
                haploid_components[accessible_pos[pos]][j] = accessible_pos[hap_cuts[j][i]]
                haploid_components[accessible_pos[pos] + 1][j] = accessible_pos[hap_cuts[j][i]]

    superreads = ReadSet()
    for i in range(ploidy):
        read = Read("superread {}".format(i + 1), 0, 0)
        # insert alleles
        for j, allele in enumerate(haplotypes[i]):
            if allele < 0:
                continue
            read.add_variant(accessible_pos[j], allele, 0)
        superreads.add(read)

    # Plot option
    if param.plot_clusters or param.plot_threading:
        timers.start("create_plots")
        draw_plots(
            block_readsets,
            clustering,
            threading,
            haplotypes,
            cuts[:-1],
            phasable_variant_table,
            param.plot_clusters,
            param.plot_threading,
            output,
        )
        timers.stop("create_plots")

    # Return results
    return components, haploid_components, superreads


def phase_single_block(block_id, block_readset, genotype_slice, param, timers):
    """
    Takes as input data the reads from a single (pre-computed) block and the genotypes for all
    variants inside the block. Runs a three-phase algorithm to compute a phasing for this isolated
    block. Input are four objects:

    block_readset -- Input reads
    genotype_slice -- Genotype (as dictionary) for every position
    param -- Object containing phasing parameters
    timers -- Timer object to measure time
    """

    block_num_vars = len(block_readset.get_positions())

    # Check for singleton blocks and handle them differently (for efficiency reasons)
    if block_num_vars == 1:

        # construct trivial solution for singleton blocks, by using the genotype as phasing
        allele_to_id = dict()
        for allele in genotype_slice[0]:
            if allele not in allele_to_id:
                allele_to_id[allele] = len(allele_to_id)
        clusts = [[] for _ in range(len(allele_to_id))]
        for i, read in enumerate(block_readset):
            clusts[allele_to_id[read[0].allele]].append(i)

        paths = [[]]
        haps = []
        for allele in genotype_slice[0]:
            for i in range(genotype_slice[0][allele]):
                paths[0].append(allele_to_id[allele])
                haps.append([allele])

        return SingleBlockResult(
            block_id, clusts, paths, [0], [[0] for _ in range(param.ploidy)], haps
        )

    # Block is non-singleton here, so run the normal routine
    # Phase I: Cluster Editing

    # Compute similarity values for all read pairs
    timers.start("read_scoring")
    logger.debug("Computing similarities for read pairs ..")
    sim = scoreReadset(block_readset, param.min_overlap, param.ploidy, 0.07)

    # Run cluster editing
    logger.debug(
        "Solving cluster editing instance with {} nodes and {} edges ..".format(
            len(block_readset), len(sim)
        )
    )
    timers.stop("read_scoring")
    timers.start("clustering")
    solver = ClusterEditingSolver(sim, param.ce_bundle_edges)
    clustering = solver.run()
    del solver
    del sim

    # Add trailing isolated nodes to single-ton clusters, if missing
    nodes_in_c = sum([len(c) for c in clustering])
    for i in range(nodes_in_c, len(block_readset)):
        clustering.append([i])

    timers.stop("clustering")

    # Phase II: Threading

    # Assemble clusters to haplotypes
    logger.debug("Threading haplotypes through {} clusters ..\r".format(len(clustering)))
    timers.start("threading")

    # Add dynamic programming for finding the most likely subset of clusters
    paths, haplotypes = run_threading(
        block_readset,
        clustering,
        param.ploidy,
        genotypes=None if param.distrust_genotypes else genotype_slice,
    )
    timers.stop("threading")

    # Phase III: Reordering

    logger.debug("Reordering ambiguous sites ..\r")
    timers.start("reordering")
    cut_positions, haploid_cuts, path, haplotypes = run_reordering(
        block_readset, clustering, paths, haplotypes, param.block_cut_sensitivity
    )
    timers.stop("reordering")

    # collect results from threading
    return SingleBlockResult(
        block_id=block_id,
        clustering=clustering,
        paths=paths,
        cuts=cut_positions,
        hap_cuts=haploid_cuts,
        haplotypes=haplotypes,
    )


def phase_single_block_mt(
    block_id, block_readset, genotype_slice, param, timers, job_id, num_blocks
):
    """
    Wrapper for the phase_single_block() function. Carries a block_id through to the results
    and unpythonizes the given readset, because cython objects are very troublesome to hardcopy.
    """
    block_vars = len(block_readset.get_positions())
    if block_vars > 1:
        logger.info(
            "Phasing block {} of {} with {} reads and {} variants.".format(
                job_id + 1, num_blocks, len(block_readset), block_vars
            )
        )
    result = phase_single_block(block_id, block_readset, genotype_slice, param, timers)
    del block_readset
    if block_vars > 1:
        logger.info("Finished block {}.".format(job_id + 1))
    return result


def aggregate_results(results, ploidy):
    """
    Collects all blockwise phasing results and aggregates them into one list for each type of
    information. Local ids and indices are converted to globals ones in this step.
    """

    clustering, cuts = [], []
    paths = []
    hap_cuts = [[] for _ in range(ploidy)]
    haplotypes = [[] for _ in range(ploidy)]
    rid_offset, cid_offset, pos_offset = 0, 0, 0
    for r in results:
        clustering += [[rid_offset + rid for rid in clust] for clust in r.clustering]
        paths += [[cid_offset + cid for cid in p] for p in r.paths]
        cuts += [pos_offset + c for c in r.cuts]
        for hap_cut, ext in zip(hap_cuts, r.hap_cuts):
            hap_cut += [pos_offset + h for h in ext]
        for hap, ext in zip(haplotypes, r.haplotypes):
            hap += ext
        rid_offset = max([rid for clust in clustering for rid in clust])
        cid_offset = len(clustering)
        pos_offset = len(haplotypes[0])

    return clustering, paths, haplotypes, cuts, hap_cuts


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
    arg(
        "--indels",
        dest="indels",
        default=False,
        action="store_true",
        help="Also phase indels (default: do not phase indels)",
    )
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


def main(args):
    run_polyphase(**vars(args))
