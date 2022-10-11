"""
Algorithmic core of whatshap polyphase

TODO: More description

"""
import logging

from itertools import chain
from multiprocessing import Pool

from whatshap.polyphase import PolyphaseBlockResult, compute_block_starts
from whatshap.polyphase.reorder import run_reordering
from whatshap.polyphase_solver import ClusterEditingSolver, scoreReadset
from whatshap.polyphase.threading import run_threading

__author__ = "Sven Schrinner"

logger = logging.getLogger(__name__)


def solve_polyphase_instance(allele_matrix, genotype_list, param, timers, quiet=False):

    # Precompute block borders based on read coverage and linkage between variants
    if not quiet:
        logger.info("Detecting connected components with weak interconnect ..")
    timers.start("detecting_blocks")
    num_vars = len(allele_matrix.getPositions())
    ploidy = param.ploidy
    if param.block_cut_sensitivity == 0:
        block_starts = [0]
    elif param.block_cut_sensitivity == 1:
        block_starts = compute_block_starts(allele_matrix, ploidy, single_linkage=True)
    else:
        block_starts = compute_block_starts(allele_matrix, ploidy, single_linkage=False)

    # Set block borders and split readset
    block_starts.append(num_vars)
    num_blocks = sum(1 for i, j in zip(block_starts[:-1], block_starts[1:]) if j > i + 1)
    if not quiet:
        logger.info(
            "Split heterozygous variants into {} blocks (and {} singleton blocks).".format(
                num_blocks, len(block_starts) - num_blocks - 1
            )
        )

    submatrices = [
        allele_matrix.extractInterval(s, e) for s, e in zip(block_starts[:-1], block_starts[1:])
    ]
    timers.stop("detecting_blocks")

    # Process blocks independently
    results = []

    # Create genotype slices for blocks
    gt_slices = []
    for block_id, submatrix in enumerate(submatrices):
        block_start = block_starts[block_id]
        block_end = block_starts[block_id + 1]
        block_num_vars = block_end - block_start

        assert submatrix.getNumPositions() == block_num_vars
        gt_slices.append(genotype_list[block_start:block_end])

    processed_blocks = 0

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
        for block_id, submatrix in enumerate(submatrices):
            block_num_vars = block_starts[block_id + 1] - block_starts[block_id]
            if block_num_vars > 1:
                processed_blocks += 1
                if not quiet:
                    logger.info(
                        "Processing block {} of {} with {} reads and {} variants.".format(
                            processed_blocks,
                            num_blocks,
                            len(submatrix),
                            block_num_vars,
                        )
                    )
            results.append(
                phase_single_block(block_id, submatrix, gt_slices[block_id], param, timers)
            )

    else:
        # sort block by descending size (4/3-approximation for scheduling problem)
        joblist = [(i, sub) for i, sub in enumerate(submatrices)]
        joblist.sort(key=lambda x: -len(x[1]))

        timers.start("phase_blocks")
        with Pool(processes=param.threads) as pool:
            process_results = [
                pool.apply_async(
                    phase_single_block_mt,
                    (
                        block_id,
                        submatrices[block_id],
                        gt_slices[block_id],
                        param,
                        timers,
                        job_id,
                        num_blocks,
                    ),
                )
                for job_id, (block_id, submatrix) in enumerate(joblist)
            ]
            # collect all blockwise results
            blockwise_results = [res.get() for res in process_results]
            results = sorted(blockwise_results, key=lambda x: x.block_id)

        timers.stop("phase_blocks")

    # Free memory
    for sub in submatrices:
        del sub

    # Aggregate blockwise results
    return aggregate_results(results, ploidy)


def phase_single_block(block_id, allele_matrix, genotype_slice, param, timers, quiet=False):
    """
    Takes as input data the reads from a single (pre-computed) block and the genotypes for all
    variants inside the block. Runs a three-phase algorithm to compute a phasing for this isolated
    block. Input are four objects:

    allele_matrix -- Input reads
    genotype_slice -- Genotype (as dictionary) for every position
    param -- Object containing phasing parameters
    timers -- Timer object to measure time
    """

    block_num_vars = allele_matrix.getNumPositions()

    # Check for singleton blocks and handle them differently (for efficiency reasons)
    if block_num_vars == 1:

        # construct trivial solution for singleton blocks, by using the genotype as phasing
        g = genotype_slice[0]
        clusts = [[i for i, r in enumerate(allele_matrix) if r and r[0][1] == a] for a in g]
        paths = [sorted(list(chain(*[[i] * g[a] for i, a in enumerate(g)])))]
        haps = sorted(list(chain(*[[[a]] * g[a] for a in g])))
        return PolyphaseBlockResult(block_id, clusts, paths, [0], [[0]] * param.ploidy, haps)

    # Block is non-singleton here, so run the normal routine
    # Phase I: Cluster Editing

    # Compute similarity values for all read pairs
    timers.start("read_scoring")
    logger.debug("Computing similarities for read pairs ..")
    sim = scoreReadset(allele_matrix, param.min_overlap, param.ploidy, 0.07)
    timers.stop("read_scoring")

    # Run cluster editing
    timers.start("clustering")
    logger.debug(
        "Solving cluster editing instance with {} nodes and {} edges ..".format(
            len(allele_matrix), len(sim)
        )
    )
    solver = ClusterEditingSolver(sim, param.ce_bundle_edges)
    clustering = solver.run()
    del solver
    del sim

    # Add trailing isolated nodes to single-ton clusters, if missing
    nodes_in_c = sum([len(c) for c in clustering])
    for i in range(nodes_in_c, len(allele_matrix)):
        clustering.append([i])

    timers.stop("clustering")

    # Phase II: Threading

    # Assemble clusters to haplotypes
    logger.debug("Threading haplotypes through {} clusters ..\r".format(len(clustering)))
    timers.start("threading")

    # Add dynamic programming for finding the most likely subset of clusters
    paths, haplotypes = run_threading(
        allele_matrix,
        clustering,
        param.ploidy,
        genotypes=None if param.distrust_genotypes else genotype_slice,
    )
    timers.stop("threading")

    # Phase III: Reordering

    logger.debug("Reordering ambiguous sites ..\r")
    timers.start("reordering")
    cut_positions, haploid_cuts, path, haplotypes = run_reordering(
        allele_matrix, clustering, paths, haplotypes, param.block_cut_sensitivity
    )
    timers.stop("reordering")

    # collect results from threading
    return PolyphaseBlockResult(
        block_id=block_id,
        clustering=clustering,
        paths=paths,
        cuts=cut_positions,
        hap_cuts=haploid_cuts,
        haplotypes=haplotypes,
    )


def phase_single_block_mt(
    block_id, allele_matrix, genotype_slice, param, timers, job_id, num_blocks, quiet=False
):
    """
    Wrapper for the phase_single_block() function. Carries a block_id through to the results
    and unpythonizes the given readset, because cython objects are very troublesome to hardcopy.
    """
    block_vars = allele_matrix.getNumPositions()
    if block_vars > 1 and not quiet:
        logger.info(
            "Phasing block {} of {} with {} reads and {} variants.".format(
                job_id + 1, num_blocks, len(allele_matrix), block_vars
            )
        )
    result = phase_single_block(block_id, allele_matrix, genotype_slice, param, timers)
    if block_vars > 1 and not quiet:
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
