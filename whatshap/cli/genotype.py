"""
Genotype variants

Runs only the genotyping algorithm. Genotype Likelihoods are computed using the
forward backward algorithm.
"""
import logging
import sys
import platform
from argparse import SUPPRESS
from typing import Sequence, Optional

from contextlib import ExitStack
from whatshap import __version__
from whatshap.vcf import VcfReader, GenotypeVcfWriter
from whatshap.core import (
    ReadSet,
    Pedigree,
    NumericSampleIds,
    PhredGenotypeLikelihoods,
    GenotypeDPTable,
    compute_genotypes,
    Genotype,
)
from whatshap.pedigree import (
    PedReader,
    UniformRecombinationCostComputer,
    GeneticMapRecombinationCostComputer,
)
from whatshap.timer import StageTimer
from whatshap.cli import log_memory_usage
from whatshap.cli.phase import select_reads, setup_families
from whatshap.cli import CommandLineError, PhasedInputReader


logger = logging.getLogger(__name__)


def int_to_diploid_biallelic_gt(numeric_repr):
    """Converts the classic numeric representation of biallelic, diploid genotypes
    into a genotype object
    """
    if numeric_repr == 0:
        return Genotype([0, 0])
    elif numeric_repr == 1:
        return Genotype([0, 1])
    elif numeric_repr == 2:
        return Genotype([1, 1])
    else:
        return Genotype([])


def determine_genotype(likelihoods: Sequence[float], threshold_prob: float) -> float:
    """given genotype likelihoods for 0/0, 0/1, 1/1, determines likeliest genotype"""

    to_sort = [
        (likelihoods[int_to_diploid_biallelic_gt(0)], 0),
        (likelihoods[int_to_diploid_biallelic_gt(1)], 1),
        (likelihoods[int_to_diploid_biallelic_gt(2)], 2),
    ]
    to_sort.sort(key=lambda x: x[0])

    # make sure there is a unique maximum which is greater than the threshold
    if (to_sort[2][0] > to_sort[1][0]) and (to_sort[2][0] > threshold_prob):
        return int_to_diploid_biallelic_gt(to_sort[2][1])
    else:
        return int_to_diploid_biallelic_gt(-1)


def run_genotype(
    phase_input_files,
    variant_file,
    reference=None,
    output=sys.stdout,
    samples=None,
    chromosomes=None,
    ignore_read_groups=False,
    only_snvs=False,
    mapping_quality=20,
    max_coverage=15,
    nopriors=False,
    ped=None,
    recombrate=1.26,
    genmap=None,
    gt_qual_threshold=0,
    prioroutput=None,
    constant=0.0,
    overhang=10,
    affine_gap=False,
    gap_start=10,
    gap_extend=7,
    mismatch=15,
    write_command_line_header=True,
    use_ped_samples=False,
    use_kmerald=False,
    kmeralign_costs_path=False,
    kmer_size=7,
    kmerald_gappenalty=40,
    kmerald_window=25,
):
    """
    For now: this function only runs the genotyping algorithm. Genotype likelihoods for
    all variants are computed using the forward backward algorithm
    """
    timers = StageTimer()
    logger.info(
        "This is WhatsHap (genotyping) %s running under Python %s",
        __version__,
        platform.python_version(),
    )
    if write_command_line_header:
        command_line = "(whatshap {}) {}".format(__version__, " ".join(sys.argv[1:]))
    else:
        command_line = None
    with ExitStack() as stack:
        # read the given input files (BAMs, VCFs, ref...)
        numeric_sample_ids = NumericSampleIds()
        phased_input_reader = stack.enter_context(
            PhasedInputReader(
                phase_input_files,
                reference,
                numeric_sample_ids,
                ignore_read_groups,
                only_snvs=only_snvs,
                mapq_threshold=mapping_quality,
                overhang=overhang,
                affine=affine_gap,
                gap_start=gap_start,
                gap_extend=gap_extend,
                default_mismatch=mismatch,
                use_kmerald=use_kmerald,
                kmeralign_costs_path=kmeralign_costs_path,
                kmer_size=kmer_size,
                kmerald_gappenalty=kmerald_gappenalty,
                kmerald_window=kmerald_window,
            )
        )
        show_phase_vcfs = phased_input_reader.has_vcfs

        # vcf writer for final genotype likelihoods
        vcf_writer = stack.enter_context(
            GenotypeVcfWriter(command_line=command_line, in_path=variant_file, out_file=output)
        )
        # vcf writer for only the prior likelihoods (if output is desired)
        prior_vcf_writer: Optional[GenotypeVcfWriter] = None
        if prioroutput is not None:
            prior_vcf_writer = stack.enter_context(
                GenotypeVcfWriter(
                    command_line=command_line,
                    in_path=variant_file,
                    out_file=stack.enter_context(open(prioroutput, "w")),
                )
            )

        # parse vcf with input variants
        # remove all likelihoods that may already be present
        vcf_reader = stack.enter_context(
            VcfReader(
                variant_file, only_snvs=only_snvs, genotype_likelihoods=False, ignore_genotypes=True
            )
        )

        if ignore_read_groups and not samples and len(vcf_reader.samples) > 1:
            raise CommandLineError(
                "When using --ignore-read-groups on a VCF with "
                "multiple samples, --sample must also be used."
            )
        if not samples:
            samples = vcf_reader.samples

        # if --use-ped-samples is set, use only samples from PED file
        if ped and use_ped_samples:
            samples = set()
            for trio in PedReader(ped):
                if trio.child is None or trio.mother is None or trio.father is None:
                    continue
                samples.add(trio.mother)
                samples.add(trio.father)
                samples.add(trio.child)

        vcf_sample_set = set(vcf_reader.samples)
        for sample in samples:
            if sample not in vcf_sample_set:
                raise CommandLineError(
                    f"Sample {sample!r} requested on command-line not found in VCF"
                )

        if ped and genmap:
            logger.info("Using region-specific recombination rates from genetic map %s.", genmap)
            recombination_cost_computer = GeneticMapRecombinationCostComputer(genmap)
        else:
            if ped:
                logger.info("Using uniform recombination rate of %g cM/Mb.", recombrate)
            recombination_cost_computer = UniformRecombinationCostComputer(recombrate)

        samples = frozenset(samples)
        families, family_trios = setup_families(samples, ped, max_coverage)
        for trios in family_trios.values():
            for trio in trios:
                # Ensure that all mentioned individuals have a numeric id
                _ = numeric_sample_ids[trio.child]

        # Read phase information provided as VCF files, if provided.
        with timers("parse_phasing_vcfs"):
            phased_input_reader.read_vcfs()

        # compute genotype likelihood threshold
        gt_prob = 1.0 - (10 ** (-gt_qual_threshold / 10.0))

        for variant_table in timers.iterate("parse_vcf", vcf_reader):
            # create a mapping of genome positions to indices
            var_to_pos = dict()
            for i in range(len(variant_table.variants)):
                var_to_pos[variant_table.variants[i].position] = i

            chromosome = variant_table.chromosome
            if (not chromosomes) or (chromosome in chromosomes):
                logger.info("======== Working on chromosome %r", chromosome)
            else:
                logger.info(
                    "Leaving chromosome %r unchanged (present in VCF but not requested by option --chromosome)",
                    chromosome,
                )
                vcf_writer.write_unchanged(chromosome)
                if prioroutput is not None:
                    prior_vcf_writer.write_unchanged(chromosome)
                continue

            positions = [v.position for v in variant_table.variants]
            if not nopriors:
                # compute prior genotype likelihoods based on all reads
                for sample in samples:
                    logger.info("---- Initial genotyping of %s", sample)
                    with timers("read_bam"):
                        readset, vcf_source_ids = phased_input_reader.read(
                            chromosome, variant_table.variants, sample, read_vcf=False
                        )
                        readset.sort()
                        genotypes, genotype_likelihoods = compute_genotypes(readset, positions)
                        # recompute genotypes based on given threshold
                        reg_genotype_likelihoods = []
                        for gl in range(len(genotype_likelihoods)):
                            norm_sum = (
                                genotype_likelihoods[gl][0]
                                + genotype_likelihoods[gl][1]
                                + genotype_likelihoods[gl][2]
                                + 3 * constant
                            )
                            regularized = PhredGenotypeLikelihoods(
                                [
                                    (genotype_likelihoods[gl][0] + constant) / norm_sum,
                                    (genotype_likelihoods[gl][1] + constant) / norm_sum,
                                    (genotype_likelihoods[gl][2] + constant) / norm_sum,
                                ]
                            )
                            genotypes[gl] = determine_genotype(regularized, gt_prob)
                            assert isinstance(genotypes[gl], Genotype)
                            reg_genotype_likelihoods.append(regularized)
                        variant_table.set_genotype_likelihoods_of(
                            sample,
                            [PhredGenotypeLikelihoods(list(gl)) for gl in reg_genotype_likelihoods],
                        )
                        variant_table.set_genotypes_of(sample, genotypes)
            else:
                # use uniform genotype likelihoods for all individuals
                for sample in samples:
                    variant_table.set_genotype_likelihoods_of(
                        sample, [PhredGenotypeLikelihoods([1 / 3, 1 / 3, 1 / 3])] * len(positions)
                    )

            # if desired, output the priors in separate vcf
            if prioroutput is not None:
                prior_vcf_writer.write_genotypes(chromosome, variant_table, only_snvs)

            # Iterate over all families to process, i.e. a separate DP table is created
            # for each family.
            for representative_sample, family in sorted(families.items()):
                if len(family) == 1:
                    logger.info("---- Processing individual %s", representative_sample)
                else:
                    logger.info("---- Processing family with individuals: %s", ",".join(family))
                max_coverage_per_sample = max(1, max_coverage // len(family))
                logger.info("Using maximum coverage per sample of %dX", max_coverage_per_sample)
                trios = family_trios[representative_sample]
                assert (len(family) == 1) or (len(trios) > 0)

                # Get the reads belonging to each sample
                readsets = dict()
                for sample in family:
                    with timers("read_bam"):
                        readset, vcf_source_ids = phased_input_reader.read(
                            chromosome, variant_table.variants, sample
                        )

                    with timers("select"):
                        readset = readset.subset(
                            [i for i, read in enumerate(readset) if len(read) >= 2]
                        )
                        logger.info(
                            "Kept %d reads that cover at least two variants each", len(readset)
                        )
                        selected_reads = select_reads(
                            readset, max_coverage_per_sample, preferred_source_ids=vcf_source_ids
                        )
                    readsets[sample] = selected_reads

                # Merge reads into one ReadSet (note that each Read object
                # knows the sample it originated from).
                all_reads = ReadSet()
                for sample, readset in readsets.items():
                    for read in readset:
                        assert read.is_sorted(), "Add a read.sort() here"
                        all_reads.add(read)

                all_reads.sort()

                # Determine which variants can (in principle) be phased
                accessible_positions = sorted(all_reads.get_positions())
                logger.info(
                    "Variants covered by at least one phase-informative "
                    "read in at least one individual after read selection: %d",
                    len(accessible_positions),
                )

                # Create Pedigree
                pedigree = Pedigree(numeric_sample_ids)
                for sample in family:
                    # genotypes are assumed to be unknown, so ignore information that
                    # might already be present in the input vcf
                    all_genotype_likelihoods = variant_table.genotype_likelihoods_of(sample)
                    genotype_l = [
                        all_genotype_likelihoods[var_to_pos[a_p]] for a_p in accessible_positions
                    ]
                    pedigree.add_individual(
                        sample, [Genotype([]) for i in range(len(accessible_positions))], genotype_l
                    )
                for trio in trios:
                    pedigree.add_relationship(
                        father_id=trio.father, mother_id=trio.mother, child_id=trio.child
                    )

                recombination_costs = recombination_cost_computer.compute(accessible_positions)

                # Finally, run genotyping algorithm
                with timers("genotyping"):
                    problem_name = "genotyping"
                    logger.info(
                        "Genotype %d sample%s by solving the %s problem ...",
                        len(family),
                        "s" if len(family) > 1 else "",
                        problem_name,
                    )
                    forward_backward_table = GenotypeDPTable(
                        numeric_sample_ids,
                        all_reads,
                        recombination_costs,
                        pedigree,
                        accessible_positions,
                    )
                    # store results
                    for s in family:
                        likelihood_list = variant_table.genotype_likelihoods_of(s)
                        genotypes_list = variant_table.genotypes_of(s)

                        for pos in range(len(accessible_positions)):
                            likelihoods = forward_backward_table.get_genotype_likelihoods(s, pos)

                            # compute genotypes from likelihoods and store information
                            geno = determine_genotype(likelihoods, gt_prob)
                            assert isinstance(geno, Genotype)
                            genotypes_list[var_to_pos[accessible_positions[pos]]] = geno
                            likelihood_list[var_to_pos[accessible_positions[pos]]] = likelihoods

                        variant_table.set_genotypes_of(s, genotypes_list)
                        variant_table.set_genotype_likelihoods_of(s, likelihood_list)

            with timers("write_vcf"):
                logger.info("======== Writing VCF")
                vcf_writer.write_genotypes(chromosome, variant_table, only_snvs)
                logger.info("Done writing VCF")

            logger.debug("Chromosome %r finished", chromosome)

    logger.info("\n== SUMMARY ==")
    total_time = timers.total()
    log_memory_usage()
    logger.info("Time spent reading BAM:                      %6.1f s", timers.elapsed("read_bam"))
    logger.info("Time spent parsing VCF:                      %6.1f s", timers.elapsed("parse_vcf"))
    if show_phase_vcfs:
        logger.info(
            "Time spent parsing input phasings from VCFs: %6.1f s",
            timers.elapsed("parse_phasing_vcfs"),
        )
    logger.info("Time spent selecting reads:                  %6.1f s", timers.elapsed("select"))
    logger.info(
        "Time spent genotyping:                          %6.1f s", timers.elapsed("genotyping")
    )
    logger.info("Time spent writing VCF:                      %6.1f s", timers.elapsed("write_vcf"))
    logger.info("Time spent on rest:                          %6.1f s", total_time - timers.sum())
    logger.info("Total elapsed time:                          %6.1f s", total_time)


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    # Positional arguments
    arg('variant_file', metavar='VCF', help='VCF file with variants to be genotyped (can be gzip-compressed)')
    arg('phase_input_files', nargs='*', metavar='PHASEINPUT',
        help='BAM or VCF file(s) with phase information, either through sequencing reads (BAM) or through phased blocks (VCF)')

    arg('-o', '--output', default=sys.stdout,
        help='Output VCF file. Add .gz to the file name to get compressed output. '
        'If omitted, use standard output.')
    arg('--reference', '-r', metavar='FASTA',
        help='Reference file. Provide this to detect alleles through re-alignment. '
        'If no index (.fai) exists, it will be created')

    arg = parser.add_argument_group('Input pre-processing, selection and filtering').add_argument
    arg('--max-coverage', '-H', metavar='MAXCOV', default=15, type=int,
        help='Reduce coverage to at most MAXCOV (default: %(default)s).')
    arg('--mapping-quality', '--mapq', metavar='QUAL',
        default=20, type=int, help='Minimum mapping quality (default: %(default)s)')
    arg('--indels', dest='indels_used', action='store_true',
        help=SUPPRESS)
    arg('--only-snvs', action='store_true', help='Genotype only SNVs')
    arg('--ignore-read-groups', default=False, action='store_true',
        help='Ignore read groups in BAM header and assume all reads come '
        'from the same sample.')
    arg('--sample', dest='samples', metavar='SAMPLE', default=[], action='append',
        help='Name of a sample to genotype. If not given, all samples in the '
        'input VCF are genotyped. Can be used multiple times.')
    arg('--chromosome', dest='chromosomes', metavar='CHROMOSOME', default=[], action='append',
        help='Name of chromosome to genotyped. If not given, all chromosomes in the '
        'input VCF are genotyped. Can be used multiple times.')
    arg('--gt-qual-threshold', metavar='GTQUALTHRESHOLD', type=float, default=0,
        help='Phred scaled error probability threshold used for genotyping (default: %(default)s). Must be at least 0. '
        'If error probability of genotype is higher, genotype ./. is output.')
    arg('--no-priors', dest='nopriors', default=False, action='store_true',
        help='Skip initial prior genotyping and use uniform priors (default: %(default)s).')
    arg('-p', '--prioroutput', default=None,
        help='output prior genotype likelihoods to the given file.')
    arg('--overhang', metavar='OVERHANG', default=10, type=int,
        help='When --reference is used, extend alignment by this many bases to left and right when realigning (default: %(default)s).')
    arg('--constant', metavar='CONSTANT', default=0, type=float,
        help='This constant is used to regularize the priors (default: %(default)s).')
    arg('--affine-gap', default=False, action='store_true',
        help='When detecting alleles through re-alignment, use affine gap costs (EXPERIMENTAL).')
    arg('--gap-start', metavar='GAPSTART', default=10, type=float,
        help='gap starting penalty in case affine gap costs are used (default: %(default)s).')
    arg('--gap-extend', metavar='GAPEXTEND', default=7, type=float,
        help='gap extend penalty in case affine gap costs are used (default: %(default)s).')
    arg('--mismatch', metavar='MISMATCH', default=15, type=float,
        help='mismatch cost in case affine gap costs are used (default: %(default)s)')

    arg = parser.add_argument_group('Pedigree genotyping').add_argument
    arg('--ped', metavar='PED/FAM',
        help='Use pedigree information in PED file to improve genotyping '
        '(switches to PedMEC algorithm). Columns 2, 3, 4 must refer to child, '
        'mother, and father sample names as used in the VCF and BAM. Other '
        'columns are ignored (EXPERIMENTAL).')
    arg('--recombrate', metavar='RECOMBRATE', type=float, default=1.26,
        help='Recombination rate in cM/Mb (used with --ped). If given, a constant recombination '
        'rate is assumed (default: %(default)gcM/Mb).')
    arg('--genmap', metavar='FILE',
        help='File with genetic map (used with --ped) to be used instead of constant recombination '
        'rate, i.e. overrides option --recombrate.')
    arg('--use-ped-samples', dest='use_ped_samples',
        action='store_true', default=False,
        help='Only work on samples mentioned in the provided PED file.')

    arg = parser.add_argument_group('kmerald based genotyping').add_argument
    arg('--use-kmerald', default=False, action='store_true',
        help='Use kmerald for detecting alleles through re-alignment.')
    arg('--kmeralign-costs', dest='kmeralign_costs_path', metavar='COSTS', default=None,
        help='Error model based costs used by kmerald during re-alignment.')
    arg('--kmer-size', metavar='KMER', type=int, default=None,
        help='kmer size used by kmerald during re-alignment.')
    arg('--kmerald-gappenalty', metavar='GAP', type=float, default=None,
        help='Gap penalty used by kmerald during re-alignment.')
    arg('--kmerald-window', metavar='WINDOW', type=float, default=None,
        help='Consider this many bases on the left and side of a variant position for kmerald based re-alignment.')
# fmt: on


def validate(args, parser):
    if args.ignore_read_groups and args.ped:
        parser.error("Option --ignore-read-groups cannot be used together with --ped")
    if args.genmap and not args.ped:
        parser.error("Option --genmap can only be used together with --ped")
    if args.genmap and len(args.chromosomes) != 1:
        parser.error(
            "Option --genmap can only be used when working on exactly one chromosome (use --chromosome)"
        )
    if len(args.phase_input_files) == 0:
        parser.error("Not providing any PHASEINPUT files not allowed for genotyping.")
    if args.gt_qual_threshold < 0:
        parser.error("Genotype quality threshold (gt-qual-threshold) must be at least 0.")
    if args.prioroutput is not None and args.nopriors:
        parser.error("Genotype priors are only computed if --no-priors is NOT set.")
    if args.constant != 0 and args.nopriors:
        parser.error("--constant can only be used if --no-priors is NOT set..")
    if args.affine_gap and not args.reference:
        parser.error("Option --affine-gap can only be used together with --reference.")
    if args.use_ped_samples and not args.ped:
        parser.error("Option --use-ped-samples can only be used when PED file is provided (--ped).")
    if args.use_ped_samples and args.samples:
        parser.error("Option --use-ped-samples cannot be used together with --samples")
    if args.indels_used:
        logger.warning("Ignoring --indels as indel genotyping is default in WhatsHap 2.0+")
    if args.use_kmerald and not args.kmeralign_costs_path:
        parser.error(
            "Option --use-kmerald can only be used when the costs to be used for kmer alignment --kmeralign-costs are provided."
        )


def main(args):
    del args.indels_used
    run_genotype(**vars(args))
