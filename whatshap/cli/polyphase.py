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
from scipy.stats import binom_test
from queue import Queue
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
    scoreReadsetLocal,
)
from whatshap.cli import log_memory_usage, PhasedInputReader, CommandLineError
from whatshap.polyphaseplots import draw_plots
from whatshap.threading import (
    run_threading,
    get_local_cluster_consensus_withfrac,
    get_position_map,
    get_pos_to_clusters_map,
    get_cluster_start_end_positions,
    get_coverage,
    get_coverage_absolute,
)
from whatshap.timer import StageTimer
from whatshap.vcf import VcfReader, PhasedVcfWriter, PloidyError

__author__ = "Jana Ebler, Sven Schrinner"

PhasingParameter = namedtuple(
    "PhasingParameter",
    [
        "ploidy",
        "verify_genotypes",
        "ce_bundle_edges",
        "min_overlap",
        "ce_refinements",
        "block_cut_sensitivity",
        "plot_clusters",
        "plot_threading",
        "threads",
    ],
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
    write_command_line_header=True,
    read_list_filename=None,
    ce_bundle_edges=False,
    min_overlap=2,
    plot_clusters=False,
    plot_threading=False,
    ce_refinements=5,
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

        if block_cut_sensitivity < 0:
            logger.warning(
                "Block cut sensitivity was set to negative value. Lowest value (0) is assumed instead."
            )
            block_cut_sensitivity = 0
        elif block_cut_sensitivity > 5:
            logger.warning(
                "Block cut sensitivity level too large. Assuming highest valid value (5) instead."
            )
            block_cut_sensitivity = 5

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
            min_overlap=min_overlap,
            ce_refinements=ce_refinements,
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
                        print(computed_genotypes, len(computed_genotypes))
                        print(given_genotypes, len(given_genotypes))
                        print(len(positions))
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
                    # TODO: Use genotype information to polish results
                    # assert len(changed_genotypes) == 0
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
    logger.info("Time spent reading BAM/CRAM:                 %6.1f s", timers.elapsed("read_bam"))
    logger.info("Time spent parsing VCF:                      %6.1f s", timers.elapsed("parse_vcf"))
    if verify_genotypes:
        logger.info(
            "Time spent verifying genotypes:              %6.1f s",
            timers.elapsed("verify_genotypes"),
        )
    logger.info(
        "Time spent detecting blocks:                 %6.1f s", timers.elapsed("detecting_blocks")
    )
    if threads == 1:
        logger.info(
            "Time spent scoring reads:                    %6.1f s", timers.elapsed("read_scoring")
        )
        logger.info(
            "Time spent solving cluster editing:          %6.1f s",
            timers.elapsed("solve_clusterediting"),
        )
        logger.info(
            "Time spent threading haplotypes:             %6.1f s", timers.elapsed("threading")
        )
    else:
        """
        TODO: The runtime measurement for the different stages does not properly for multithreading,
        because the global timer is not visible from within the phase_single_block_mt method.
        Workaround is to only report the total phasing time.
        """
        logger.info(
            "Time spent phasing blocks:                   %6.1f s", timers.elapsed("phase_blocks")
        )
    if plot_clusters or plot_threading:
        logger.info(
            "Time spent creating plots:                   %6.1f s", timers.elapsed("create_plots")
        )
    logger.info("Time spent writing VCF:                      %6.1f s", timers.elapsed("write_vcf"))
    logger.info(
        "Time spent on rest:                          %6.1f s", timers.total() - timers.sum()
    )
    logger.info("Total elapsed time:                          %6.1f s", timers.total())


def phase_single_individual(readset, phasable_variant_table, sample, phasing_param, output, timers):

    # Compute the genotypes that belong to the variant table and create a list of all genotypes
    genotype_list = create_genotype_list(phasable_variant_table, sample)

    # Select reads, only for debug
    # selected_reads = select_reads(readset, 120, preferred_source_ids = vcf_source_ids)
    # readset = selected_reads

    # Precompute block borders based on read coverage and linkage between variants
    logger.info("Detecting connected components with weak interconnect ..")
    timers.start("detecting_blocks")
    index, rev_index = get_position_map(readset)
    num_vars = len(rev_index)
    if phasing_param.block_cut_sensitivity == 0:
        block_starts = [0]
    elif phasing_param.block_cut_sensitivity == 1:
        block_starts = compute_linkage_based_block_starts(
            readset, index, phasing_param.ploidy, single_linkage=True
        )
    else:
        block_starts = compute_linkage_based_block_starts(
            readset, index, phasing_param.ploidy, single_linkage=False
        )

    # Set block borders and split readset
    ext_block_starts = block_starts + [num_vars]
    num_non_singleton_blocks = len(
        [i for i in range(len(block_starts)) if ext_block_starts[i] < ext_block_starts[i + 1] - 1]
    )
    logger.info(
        "Split heterozygous variants into {} blocks (and {} singleton blocks).".format(
            num_non_singleton_blocks, len(block_starts) - num_non_singleton_blocks
        )
    )

    block_readsets = split_readset(readset, ext_block_starts, index)
    timers.stop("detecting_blocks")

    # Process blocks independently
    (
        blockwise_clustering,
        blockwise_paths,
        blockwise_haplotypes,
        blockwise_cut_positions,
        blockwise_haploid_cuts,
    ) = ([], [], [], [], [])

    # Create genotype slices for blocks
    genotype_slices = []
    for block_id, block_readset in enumerate(block_readsets):
        block_start = ext_block_starts[block_id]
        block_end = ext_block_starts[block_id + 1]
        block_num_vars = block_end - block_start

        assert len(block_readset.get_positions()) == block_num_vars
        genotype_slices.append(genotype_list[block_start:block_end])

    processed_non_singleton_blocks = 0
    # use process pool for multiple threads
    if phasing_param.threads == 1:
        # for single-threading, process everything individually to minimize memory footprint
        for block_id, block_readset in enumerate(block_readsets):
            block_num_vars = ext_block_starts[block_id + 1] - ext_block_starts[block_id]
            if block_num_vars > 1:
                # Only print for non-singleton block
                processed_non_singleton_blocks += 1
                logger.info(
                    "Processing block {} of {} with {} reads and {} variants.".format(
                        processed_non_singleton_blocks,
                        num_non_singleton_blocks,
                        len(block_readset),
                        block_num_vars,
                    )
                )

            clustering, path, haplotypes, cut_positions, haploid_cuts = phase_single_block(
                block_readset, genotype_slices[block_id], phasing_param, timers
            )

            blockwise_clustering.append(clustering)
            blockwise_paths.append(path)
            blockwise_haplotypes.append(haplotypes)
            blockwise_cut_positions.append(cut_positions)
            blockwise_haploid_cuts.append(haploid_cuts)

    else:
        # sort block readsets in descending order by number of reads
        joblist = [(i, len(block_readsets[i])) for i in range(len(block_readsets))]
        joblist.sort(key=lambda x: -x[1])

        timers.start("phase_blocks")

        # process large jobs first, 4/3-approximation for scheduling problem
        with Pool(processes=phasing_param.threads) as pool:
            """
            TODO: Python's multiprocessing makes hard copies of the passed
            arguments, which is not trivial for cython objects, especially when they
            contain pointers to other cython objects. Any passed object must be
            (de)serializable (in Python: pickle).
            All other objects created in the main thread are also accessible by the
            workers, but they are handled via the copy-on-write policy. This means,
            that e.g. the large main readset is not hardcopied for every thread,
            as long as it is not modified there. Since this would cause a massive
            waste of memory, this must not be done and the main readset must
            also never be passed as argument to the workers.
            """
            process_results = [
                pool.apply_async(
                    phase_single_block_mt,
                    (
                        block_readsets[block_id],
                        genotype_slices[block_id],
                        phasing_param,
                        timers,
                        block_id,
                        job_id,
                        num_non_singleton_blocks,
                    ),
                )
                for job_id, (block_id, block_readset) in enumerate(joblist)
            ]
            blockwise_results = [res.get() for res in process_results]

            # reorder results again
            blockwise_results.sort(key=lambda x: x[-1])

            # collect all blockwise results
            for (
                clustering,
                path,
                haplotypes,
                cut_positions,
                haploid_cuts,
                block_id,
            ) in blockwise_results:
                blockwise_clustering.append(clustering)
                blockwise_paths.append(path)
                blockwise_haplotypes.append(haplotypes)
                blockwise_cut_positions.append(cut_positions)
                blockwise_haploid_cuts.append(haploid_cuts)

        timers.stop("phase_blocks")

    # Aggregate blockwise results
    clustering, threading, haplotypes, cut_positions, haploid_cuts = aggregate_phasing_blocks(
        block_starts,
        block_readsets,
        blockwise_clustering,
        blockwise_paths,
        blockwise_haplotypes,
        blockwise_cut_positions,
        blockwise_haploid_cuts,
        phasing_param,
    )

    # Summarize data for VCF file
    accessible_positions = sorted(readset.get_positions())
    components = {}
    haploid_components = {}

    ext_cuts = cut_positions + [num_vars]
    for i, cut_pos in enumerate(cut_positions):
        for pos in range(ext_cuts[i], ext_cuts[i + 1]):
            components[accessible_positions[pos]] = accessible_positions[ext_cuts[i]]
            components[accessible_positions[pos] + 1] = accessible_positions[ext_cuts[i]]
            haploid_components[accessible_positions[pos]] = [0] * phasing_param.ploidy
            haploid_components[accessible_positions[pos] + 1] = [0] * phasing_param.ploidy

    for j in range(phasing_param.ploidy):
        ext_cuts = haploid_cuts[j] + [num_vars]
        for i, cut_pos in enumerate(haploid_cuts[j]):
            for pos in range(ext_cuts[i], ext_cuts[i + 1]):
                haploid_components[accessible_positions[pos]][j] = accessible_positions[ext_cuts[i]]
                haploid_components[accessible_positions[pos] + 1][j] = accessible_positions[
                    ext_cuts[i]
                ]

    superreads = ReadSet()
    for i in range(phasing_param.ploidy):
        read = Read("superread {}".format(i + 1), 0, 0)
        # insert alleles
        for j, allele in enumerate(haplotypes[i]):
            if allele == "n":
                continue
            allele = int(allele)
            # TODO: Needs changes for multi-allelic and we might give an actual quality value
            read.add_variant(accessible_positions[j], allele, 0)
        superreads.add(read)

    # Plot option
    if phasing_param.plot_clusters or phasing_param.plot_threading:
        timers.start("create_plots")
        draw_plots(
            block_readsets,
            clustering,
            threading,
            haplotypes,
            cut_positions,
            genotype_list,
            phasable_variant_table,
            phasing_param.plot_clusters,
            phasing_param.plot_threading,
            output,
        )
        timers.stop("create_plots")

    # Return results
    return components, haploid_components, superreads


def create_genotype_list(phasable_variant_table, sample):
    """
    Creates a list, which stores a dictionary for every position. The dictionary maps every allele to its frequency in the genotype of the respective position.
    """
    all_genotypes = phasable_variant_table.genotypes_of(sample)
    genotype_list = []
    for pos in range(len(all_genotypes)):
        allele_count = dict()
        for allele in all_genotypes[pos].as_vector():
            if allele not in allele_count:
                allele_count[allele] = 0
            allele_count[allele] += 1
        genotype_list.append(allele_count)
    return genotype_list


def split_readset(readset, ext_block_starts, index):
    """
    Creates one sub-readset for every block. Reads which cross block borders are also split, so parts of them
    appear in multiple blocks. Reads inside a sub-readset are trimmed, such that they do not contain variants
    outside of their associated blocks.
    """

    var_to_block = [0 for _ in range(ext_block_starts[-1])]
    for i in range(len(ext_block_starts) - 1):
        for var in range(ext_block_starts[i], ext_block_starts[i + 1]):
            var_to_block[var] = i

    block_readsets = [ReadSet() for i in range(len(ext_block_starts) - 1)]
    for i, read in enumerate(readset):
        if not read.is_sorted():
            read.sort()
        start = var_to_block[index[read[0].position]]
        end = var_to_block[index[read[-1].position]]
        if start == end:
            # if read lies entirely in one block, copy it into according readset
            block_readsets[start].add(read)
        else:
            # split read by creating one new read for each covered block
            current_block = start
            read_slice = Read(
                name=read.name,
                source_id=read.source_id,
                sample_id=read.sample_id,
                reference_start=read.sample_id,
                BX_tag=read.BX_tag,
            )
            for variant in read:
                if var_to_block[index[variant.position]] != current_block:
                    block_readsets[current_block].add(read_slice)
                    current_block = var_to_block[index[variant.position]]
                    read_slice = Read(
                        name=str(current_block) + "_" + read.name,
                        source_id=read.source_id,
                        sample_id=read.sample_id,
                        reference_start=read.sample_id,
                        BX_tag=read.BX_tag,
                    )
                read_slice.add_variant(variant.position, variant.allele, variant.quality)
            block_readsets[current_block].add(read_slice)
    return block_readsets


def phase_single_block(block_readset, genotype_slice, phasing_param, timers):
    """
    Takes as input data the reads from a single (pre-computed) block and the genotypes for all variants inside the block.
    Also requires a ploidy and block cut sensitivity as parameters. Runs a two-phase algorithm to compute a phasing for
    this isolated block. The phasing algorithm may create additional block cut, i.e. it may indicate where new phasing
    blocks should start.

    Output are four objects:

    clustering -- A list of clusters, which were the result of the cluster editing step. A cluster is a list of read ids.
    path -- A list with one entry per variant. Each entry contains the tuple of clusters, through which the haplotypes
            have been threaded in the threading step.
    haplotypes -- A list containing the haplotypes. A haplotype is a string over allele ids, one allele per variant.
    cut_positions -- A list of variant positions, where the phasing blocks start. The positions do not refer to genome
                     positions (in base pairs), but to the variant id inside this block. 0 is always contained, since a
                     block needs to start at the beginning of the given input, but it can contain more positions.
    """

    block_num_vars = len(block_readset.get_positions())

    # Check for singleton blocks and handle them differently (for efficiency reasons)
    if block_num_vars == 1:

        # construct trivial solution for singleton blocks, by basically using the genotype as phasing
        allele_to_id = dict()
        for allele in genotype_slice[0]:
            if allele not in allele_to_id:
                allele_to_id[allele] = len(allele_to_id)
        clustering = [[] for _ in range(len(allele_to_id))]
        for i, read in enumerate(block_readset):
            clustering[allele_to_id[read[0].allele]].append(i)

        path = [[]]
        haplotypes = []
        for allele in genotype_slice[0]:
            for i in range(genotype_slice[0][allele]):
                path[0].append(allele_to_id[allele])
                haplotypes.append(str(allele))

        return clustering, path, haplotypes, [0], [[0] for _ in range(phasing_param.ploidy)]

    # Block is non-singleton here, so run the normal routine
    # Phase I: Cluster Editing

    # Compute similarity values for all read pairs
    timers.start("read_scoring")
    logger.debug("Computing similarities for read pairs ...")
    similarities = scoreReadsetLocal(block_readset, phasing_param.min_overlap, phasing_param.ploidy)

    # Run cluster editing
    logger.debug(
        "Solving cluster editing instance with {} nodes and {} edges ..".format(
            len(block_readset), len(similarities)
        )
    )
    timers.stop("read_scoring")
    timers.start("solve_clusterediting")
    solver = ClusterEditingSolver(similarities, phasing_param.ce_bundle_edges)
    clustering = solver.run()
    del solver

    # Refine clusters by solving inconsistencies in consensus
    runs_remaining = phasing_param.ce_refinements
    last_inc_count = len(clustering) * (block_num_vars)  # worst case number
    refine = True
    while refine and runs_remaining > 0:
        """
        Inconsistencies are positions, whre a cluster has a very ambiguous consensus, indicating that it contains reads from
        two or more haplotypes, which differ on some SNP variants
        """
        refine = False
        runs_remaining -= 1
        new_inc_count, seperated_reads = find_inconsistencies(
            block_readset, clustering, phasing_param.ploidy
        )
        for (r0, r1) in seperated_reads:
            similarities.set(r0, r1, -float("inf"))

        if 0 < new_inc_count < last_inc_count:
            logger.debug(
                "{} inconsistent variants found. Refining clusters ..\r".format(new_inc_count)
            )
            solver = ClusterEditingSolver(similarities, phasing_param.ce_bundle_edges)
            clustering = solver.run()
            del solver

    # Deallocate big datastructures, which are not needed anymore
    del similarities

    # Add trailing isolated nodes to single-ton clusters, if missing
    nodes_in_c = sum([len(c) for c in clustering])
    for i in range(nodes_in_c, len(block_readset)):
        clustering.append([i])

    timers.stop("solve_clusterediting")

    # Phase II: Threading

    # Assemble clusters to haplotypes
    logger.debug("Threading haplotypes through {} clusters..\r".format(len(clustering)))
    timers.start("threading")

    # Add dynamic programming for finding the most likely subset of clusters
    cut_positions, haploid_cuts, path, haplotypes = run_threading(
        block_readset,
        clustering,
        phasing_param.ploidy,
        genotype_slice,
        phasing_param.block_cut_sensitivity,
    )
    timers.stop("threading")

    # collect results from threading
    return clustering, path, haplotypes, cut_positions, haploid_cuts


def phase_single_block_mt(
    block_readset, genotype_slice, phasing_param, timers, block_id, job_id, num_blocks
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
    clustering, path, haplotypes, cut_positions, haploid_cuts = phase_single_block(
        block_readset, genotype_slice, phasing_param, timers
    )
    del block_readset
    if block_vars > 1:
        logger.info("Finished block {}.".format(job_id + 1))
    return clustering, path, haplotypes, cut_positions, haploid_cuts, block_id


def aggregate_phasing_blocks(
    block_starts,
    block_readsets,
    blockwise_clustering,
    blockwise_paths,
    blockwise_haplotypes,
    blockwise_cut_positions,
    blockwise_haploid_cuts,
    phasing_param,
):
    """
    Collects all blockwise phasing results and aggregates them into one list for each type of information. Local ids and indices are
    converted to globals ones in this step.
    """

    clustering = []
    read_id_offset = 0
    for i in range(len(block_starts)):
        for cluster in blockwise_clustering[i]:
            clustering.append(tuple([read_id + read_id_offset for read_id in cluster]))
        read_id_offset += len(block_readsets[i])

    threading = []
    c_id_offset = 0
    for i in range(len(block_starts)):
        for c_tuple in blockwise_paths[i]:
            threading.append(tuple([c_id + c_id_offset for c_id in c_tuple]))
        c_id_offset += len(blockwise_clustering[i])

    haplotypes = []
    for i in range(phasing_param.ploidy):
        haplotypes.append("".join([block[i] for block in blockwise_haplotypes]))

    cut_positions = []
    for i in range(len(block_starts)):
        for cut_pos in blockwise_cut_positions[i]:
            cut_positions.append(cut_pos + block_starts[i])

    haploid_cuts = [[] for _ in range(phasing_param.ploidy)]
    for i in range(len(block_starts)):
        for j in range(phasing_param.ploidy):
            for cut_pos in blockwise_haploid_cuts[i][j]:
                haploid_cuts[j].append(cut_pos + block_starts[i])

    return clustering, threading, haplotypes, cut_positions, haploid_cuts


def find_inconsistencies(readset, clustering, ploidy):
    # Returns the number of cluster positions with inconsistencies
    # (counts position multiple times, if multiple clusters are inconsistent there)
    # Also returns a list of read pairs, which need to be seperated
    num_inconsistent_positions = 0
    separated_pairs = []
    exp_error = 0.05
    p_val_threshold = 0.02

    # Compute consensus and coverage
    index, rev_index = get_position_map(readset)
    num_vars = len(rev_index)

    coverage = get_coverage(readset, clustering, index)
    cov_map = get_pos_to_clusters_map(coverage, ploidy)
    positions = get_cluster_start_end_positions(readset, clustering, index)
    abs_coverage = get_coverage_absolute(readset, clustering, index)
    consensus = get_local_cluster_consensus_withfrac(readset, clustering, cov_map, positions)

    # Search for positions in clusters with ambivalent consensus
    for pos in range(num_vars):
        # print(str(pos)+" -> "+str(len(coverage[pos]))+" , "+str(len(consensus[pos])))
        for c_id in coverage[pos]:
            if c_id not in consensus[pos]:
                continue
            # do binomial hypothesis test, whether the deviations from majority allele is significant enough for splitting
            abs_count = abs_coverage[pos][c_id]
            abs_deviations = int(abs_count * (1 - consensus[pos][c_id][1]))
            p_val = binom_test(abs_deviations, abs_count, exp_error, alternative="greater")
            if p_val < p_val_threshold:
                # print("   inconsistency in cluster "+str(c_id)+" at position"+str(pos)+" with coverage "+str(coverage[pos][c_id])+" and consensus "+str(consensus[pos][c_id]))
                num_inconsistent_positions += 1
                zero_reads = []
                one_reads = []
                for read in clustering[c_id]:
                    for var in readset[read]:
                        if index[var.position] == pos:
                            if var.allele == 0:
                                zero_reads.append(read)
                            else:
                                one_reads.append(read)
                for r0 in zero_reads:
                    for r1 in one_reads:
                        separated_pairs.append((r0, r1))

    return num_inconsistent_positions, separated_pairs


def compute_linkage_based_block_starts(readset, pos_index, ploidy, single_linkage=False):
    """
    Based on the connectivity of the reads, we want to divide the phasing input, as non- or poorly connected
    regions can be phased independently. This is done based on how pairs of variants are connected. There are
    two modes how to decide whether two variants are connected:

    single_linkage=True -- If there exists a read in the readset, which covers both variants, they are connected
    single_linkage=False -- In order to connect two variants, we need at least reads from ploidy-1 different
                            haplotypes. Two variants count as connected, if there sufficiently many reads covering
                            both variants, with "sufficient" meaning, that the connecting reads have a chance of
                            at least 98% that they cover at least ploidy-1 haplotypes.

    First, only consecutive pairs are inspected. Then, this connectivity is made transitive, i.e. if the pair
    (A,C) is connected, as well as the pair (B,C), then (A,B) is also connected. If the special case occurs, that
    for three variants (in this order) A and C are connected, but neither is connected to variant B in between them,
    then the variants are still divided as A|B|C, even though A and C are actually connected. This is because the
    following steps require the splits to be intervals with no "holes" inside them.
    """

    num_vars = len(pos_index)

    # special case
    if num_vars == 0:
        return []

    # cut threshold
    if ploidy == 2 or single_linkage:
        cut_threshold = 1
    else:
        cut_threshold = ploidy * ploidy
        for i in range(ploidy - 1, ploidy * ploidy):
            # chance to cover at most ploidy-2 haplotypes with i reads
            cut_threshold = i
            if ploidy * pow((ploidy - 2) / ploidy, i) < 0.02:
                cut_threshold = i
                break
    logger.debug("Cut position threshold: coverage >= {}".format(cut_threshold))

    # start by looking at neighbouring
    link_to_next = [0 for i in range(num_vars)]
    starts = []
    ends = []
    for read in readset:
        pos_list = [pos_index[var.position] for var in read]
        starts.append(pos_list[0])
        ends.append(pos_list[-1])
        for i in range(len(pos_list) - 1):
            if pos_list[i] + 1 == pos_list[i + 1]:
                link_to_next[pos_list[i]] += 1

    pos_clust = [0 for i in range(num_vars)]
    for i in range(1, num_vars):
        if link_to_next[i - 1] >= cut_threshold:
            pos_clust[i] = pos_clust[i - 1]
        else:
            pos_clust[i] = pos_clust[i - 1] + 1
    num_clust = pos_clust[-1] + 1

    # find linkage between clusters
    link_coverage = [dict() for i in range(num_clust)]
    for i, (start, end) in enumerate(zip(starts, ends)):
        covered_pos_clusts = set([pos_clust[pos_index[var.position]] for var in readset[i]])
        for p1 in covered_pos_clusts:
            for p2 in covered_pos_clusts:
                if p2 not in link_coverage[p1]:
                    link_coverage[p1][p2] = 0
                link_coverage[p1][p2] += 1

    # merge clusters
    merged_clust = [-1 for i in range(num_clust)]
    new_num_clust = 0
    for i in range(num_clust):
        if merged_clust[i] >= 0:
            continue
        q = Queue()
        q.put(i)
        merged_clust[i] = new_num_clust
        while not q.empty():
            cur = q.get()
            for linked in link_coverage[cur]:
                if merged_clust[linked] < 0 and link_coverage[cur][linked] >= cut_threshold:
                    q.put(linked)
                    merged_clust[linked] = new_num_clust
        new_num_clust += 1

    # determine cut positions
    cuts = [0]
    for i in range(1, num_vars):
        if merged_clust[pos_clust[i]] != merged_clust[pos_clust[i - 1]]:
            cuts.append(i)

    return cuts


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
    arg(
        "--output-read-list",
        metavar="FILE",
        default=None,
        dest="read_list_filename",
        help="Write reads that have been used for phasing to FILE.",
    )

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
        "--verify-genotypes",
        default=False,
        action="store_true",
        help="Verify input genotypes by re-typing them using the given reads.",
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
        "--ce-refinements",
        metavar="REFINEMENTS",
        type=int,
        default=1,
        help="Maximum number of refinement steps for internal read clustering stage (default: %(default)s).",
    )
    arg(
        "--block-cut-sensitivity",
        "-B",
        metavar="SENSITIVITY",
        type=int,
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


def validate(args, parser):
    pass


def main(args):
    run_polyphase(**vars(args))
