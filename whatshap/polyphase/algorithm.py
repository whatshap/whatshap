"""
Algorithmic core of whatshap polyphase. Split inputs into independent blocks if possible and can be
executed in parallel on a block level. Each block is processed in three phases: Read clustering,
haplotype threading and reordering.
"""
import logging

from itertools import chain
from multiprocessing import Pool
from math import log
from typing import List
from copy import copy

from whatshap.polyphase import (
    PolyphaseResult,
    PolyphaseBlockResult,
    PhaseBreakpoint,
    compute_block_starts,
)
from whatshap.polyphase.reorder import find_subinstances, integrate_sub_results, run_reordering
from whatshap.polyphase.solver import ClusterEditingSolver, scoreReadset
from whatshap.polyphase.threading import run_threading

__author__ = "Sven Schrinner"

logger = logging.getLogger(__name__)


def solve_polyphase_instance(
    allele_matrix, genotype_list, param, timers, partial_phasing=None, quiet=False
):
    """
    Entry point for polyploid phasing instances. Inputs are an allele matrix and genotypes for each
    position, among some parameters.
    """
    num_vars = len(allele_matrix.getPositions())

    assert num_vars > 0
    assert len(allele_matrix) > 0

    # Precompute block borders based on read coverage and linkage between variants
    if not quiet:
        logger.info("Detecting connected components with weak interconnect ..")
    timers.start("detecting_blocks")

    ploidy = param.ploidy
    sl = param.block_cut_sensitivity <= 1
    block_starts = compute_block_starts(allele_matrix, ploidy, single_linkage=sl)

    # Set block borders and split readset
    block_starts.append(num_vars)
    assert block_starts == sorted(list(set(block_starts)))
    num_blocks = sum(1 for i, j in zip(block_starts[:-1], block_starts[1:]) if j > i + 1)
    if not quiet:
        logger.info(
            f"Split heterozygous variants into {num_blocks} blocks (and {len(block_starts) - num_blocks - 1} singleton blocks)."
        )

    # Process blocks independently
    results = []
    processed_blocks = 0
    timers.stop("detecting_blocks")

    """
    Python's multiprocessing makes hard copies of the passed arguments, which is not trivial for
    cython objects, especially when they contain pointers to other cython objects. Any passed
    object must be (de)serializable (in Python: pickle). All other objects created in the main
    thread are also accessible by the workers, but they are handled via the copy-on-write policy.
    This means, that e.g. the large main matrix is not hardcopied for every thread, as long as it
    is not modified there. This must be ensured to prevent a massive waste of memory consumption.
    """
    if param.threads == 1:
        # for single-threading, process everything individually to minimize memory footprint
        for block_id, (start, end) in enumerate(zip(block_starts[:-1], block_starts[1:])):
            submatrix = allele_matrix.extractInterval(start, end)
            subphasing = partial_phasing.extractInterval(start, end) if partial_phasing else None
            if end - start > 1:
                processed_blocks += 1
                if not quiet:
                    logger.info(
                        f"Processing block {processed_blocks} of {num_blocks} with {len(submatrix)} reads and {end - start} variants."
                    )
            results.append(
                phase_single_block(
                    block_id, submatrix, genotype_list[start:end], subphasing, param, timers, quiet
                )
            )
            del submatrix

    else:
        # sort block by descending size (4/3-approximation for scheduling problem)
        timers.start("phase_blocks")
        joblist = list(zip(range(len(block_starts)), block_starts[:-1], block_starts[1:]))
        joblist.sort(key=lambda x: x[1] - x[2])

        with Pool(processes=param.threads) as pool:
            process_results = [
                pool.apply_async(
                    phase_single_block_mt,
                    (
                        allele_matrix,
                        partial_phasing,
                        block_id,
                        start,
                        end,
                        genotype_list[start:end],
                        param,
                        timers,
                        job_id,
                        num_blocks,
                        quiet,
                    ),
                )
                for job_id, (block_id, start, end) in enumerate(joblist)
            ]
            # collect all blockwise results
            blockwise_results = [res.get() for res in process_results]
            results = sorted(blockwise_results, key=lambda x: x.block_id)

        timers.stop("phase_blocks")

    # Aggregate blockwise results
    if partial_phasing and param.block_cut_sensitivity == 0:
        # For lowest sensitivity, do not add block starts to global breakpoint list
        # (unless the partial phasing is also interrupted there)
        borders = {partial_phasing.getFirstPos(i) for i in range(len(partial_phasing))}
    else:
        borders = []
    return aggregate_results(results, ploidy, borders)


def phase_single_block(block_id, allele_matrix, genotypes, prephasing, param, timers, quiet=False):
    """
    Takes as input data the reads from a single (pre-computed) block and the genotypes for all
    variants inside the block. Runs a three-phase algorithm to compute a phasing for this isolated
    block. Input are four objects:

    allele_matrix -- Input reads
    genotypes -- Genotype (as dictionary) for every position
    prephasing -- Positions for which a phasing already exists, encoded as reads
    param -- Object containing phasing parameters
    timers -- Timer object to measure time
    quiet -- If set, suppresses logger info output
    """

    # Check for empty/singleton blocks and handle them differently (for efficiency reasons)
    block_num_vars = allele_matrix.getNumPositions()
    if block_num_vars < 2:
        # construct trivial solution for singleton blocks, by using the genotype as phasing
        g = genotypes[0]
        clusts = [[i for i, r in enumerate(allele_matrix) if r and r[0][1] == a] for a in g]
        threads = [sorted(list(chain(*[[i] * g[a] for i, a in enumerate(g)])))]
        haps = sorted(list(chain(*[[[a]] * g[a] for a in g])))
        return PolyphaseBlockResult(block_id, clusts, threads, haps, [])

    # Block is non-singleton here, so run the normal routine
    # Phase I: Cluster Editing

    # Compute similarity values for all read pairs
    assert len(allele_matrix) > 0
    timers.start("read_scoring")
    logger.debug("Computing similarities for read pairs ..")
    sim = scoreReadset(allele_matrix, param.min_overlap, param.ploidy, 0.07)
    timers.stop("read_scoring")

    # Run cluster editing
    timers.start("clustering")
    logger.debug(
        f"Solving cluster editing instance with {len(allele_matrix)} nodes and {len(sim)} edges .."
    )
    solver = ClusterEditingSolver(sim, param.ce_bundle_edges)
    clustering = solver.run()
    del solver
    del sim

    # Add trailing isolated nodes to single-ton clusters, if missing
    nodes_in_c = sum(len(c) for c in clustering)
    for i in range(nodes_in_c, len(allele_matrix)):
        clustering.append([i])

    timers.stop("clustering")

    # Phase II: Threading

    # Assemble clusters to haplotypes
    logger.debug(f"Threading haplotypes through {len(clustering)} clusters ..\r")
    timers.start("threading")

    # Add dynamic programming for finding the most likely subset of clusters
    threads, haplotypes = run_threading(
        allele_matrix,
        clustering,
        param.ploidy,
        genotypes,
        distrust_genotypes=param.distrust_genotypes,
    )
    timers.stop("threading")

    # Phase III: Reordering

    logger.debug("Reordering ambiguous sites ..\r")
    timers.start("reordering")

    # Recursively resolve collapsed regions in clusters
    sub_instances = find_subinstances(allele_matrix, clustering, threads, haplotypes)
    sub_results = []
    sub_param = copy(param)
    sub_param.ignore_phasings = True
    sub_param.threads = 1
    for cid, thread_set, subm in sub_instances:
        assert len(subm) > 0
        snps = [allele_matrix.globalToLocal(gpos) for gpos in subm.getPositions()]
        assert all([0 <= pos < allele_matrix.getNumPositions() for pos in snps])
        subhaps = [[haplotypes[i][pos] for i in thread_set] for pos in snps]
        subgeno = [{a: h.count(a) for a in h} for h in subhaps]
        sub_param.ploidy = len(thread_set)
        timers.stop("reordering")
        res = solve_polyphase_instance(subm, subgeno, sub_param, timers, quiet=True)
        timers.start("reordering")
        sub_results.append(res)

    # collect breakpoints of sub-instances and overall instance. Update threads/haplotypes
    breakpoints = integrate_sub_results(
        allele_matrix, sub_instances, sub_results, threads, haplotypes
    )
    del sub_instances
    del sub_results

    # reorder pieces
    run_reordering(allele_matrix, clustering, threads, haplotypes, breakpoints, prephasing)

    timers.stop("reordering")

    # collect results from threading
    return PolyphaseBlockResult(
        block_id=block_id,
        clustering=[[allele_matrix.getGlobalId(r) for r in c] for c in clustering],
        threads=threads,
        haplotypes=haplotypes,
        breakpoints=breakpoints,
    )


def phase_single_block_mt(
    allele_matrix,
    partial_phasing,
    block_id,
    start,
    end,
    genotype_slice,
    param,
    timers,
    job_id,
    num_blocks,
    quiet=False,
):
    """
    Wrapper for the phase_single_block() function. Carries a block_id through to the results.
    Creates a local submatrix without modifying the given allele matrix
    """
    submatrix = allele_matrix.extractInterval(start, end)
    subphasing = partial_phasing.extractInterval(start, end) if partial_phasing else None
    block_vars = submatrix.getNumPositions()
    if block_vars > 1 and not quiet:
        logger.info(
            f"Phasing block {job_id + 1} of {num_blocks} with {len(submatrix)} reads and {block_vars} variants."
        )

    result = phase_single_block(
        block_id, submatrix, genotype_slice, subphasing, param, timers, quiet
    )
    if block_vars > 1 and not quiet:
        logger.info(f"Finished block {job_id + 1}.")
    return result


def aggregate_results(results: List[PolyphaseBlockResult], ploidy: int, borders: List[int]):
    """
    Collects all blockwise phasing results and aggregates them into one list for each type of
    information. Local ids and indices are converted to globals ones in this step.
    """
    clustering, threads, breakpoints = [], [], []
    haplotypes = [[] for _ in range(ploidy)]
    cid_offset, pos_offset = 0, 0
    for r in results:
        clustering += [clust for clust in r.clustering]
        threads += [[cid_offset + cid for cid in p] for p in r.threads]
        for hap, ext in zip(haplotypes, r.haplotypes):
            hap += ext
        # Add the start of a block as breakpoint, unless a partial phasing bridges the blocks
        if not borders or pos_offset in borders or pos_offset == 0:
            breakpoints.append(PhaseBreakpoint(pos_offset, list(range(ploidy)), 0.0))
        breakpoints += [
            PhaseBreakpoint(b.position + pos_offset, b.haplotypes, b.confidence)
            for b in r.breakpoints
        ]
        cid_offset = len(clustering)
        pos_offset = len(haplotypes[0])

    return PolyphaseResult(clustering, threads, haplotypes, breakpoints)


def compute_cut_positions(
    breakpoints: List[PhaseBreakpoint], ploidy: int, block_cut_sensitivity: int
):
    """
    Computes the cut positions for phasing blocks, based on the computed breakpoints of the
    reordering stage and the requeted block cut sensitivity.
    """

    cuts = []
    hap_cuts = [[] for _ in range(ploidy)]
    thresholds = [-float("inf"), -float("inf"), log(0.5), log(0.5), log(0.99), 0]
    thresholds_num = [ploidy, ploidy, min(ploidy, 3), 2, 2, 0]
    threshold = thresholds[block_cut_sensitivity]
    threshold_num = thresholds_num[block_cut_sensitivity]

    remaining_conf = [0.0 for _ in range(ploidy)]
    for b in breakpoints:
        # avoid duplicate cut positions
        if cuts and cuts[-1] == b.position:
            continue
        # for zero confidence, always cut
        if b.confidence == 0.0:
            cuts.append(b.position)
            for h in range(ploidy):
                hap_cuts[h].append(b.position)
            remaining_conf = [0.0 for _ in range(ploidy)]
            continue
        else:
            for h in b.haplotypes:
                remaining_conf[h] += log(b.confidence)
        if sum([1 for i in range(ploidy) if remaining_conf[i] <= threshold]) >= threshold_num:
            cuts.append(b.position)
            for h in b.haplotypes:
                hap_cuts[h].append(b.position)
            remaining_conf = [0.0 for _ in range(ploidy)]

    return cuts, hap_cuts
