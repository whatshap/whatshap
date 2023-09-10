"""
Tag reads by haplotype

Sequencing reads are read from file ALIGNMENTS (in BAM or CRAM format) and tagged reads
are written to stdout.
"""
import logging
import sys
import pysam
import hashlib
from collections import defaultdict
from typing import List, Optional, Union, Dict, Tuple, FrozenSet, Sequence, TextIO

from xopen import xopen

from contextlib import ExitStack
from whatshap import __version__
from whatshap.cli import PhasedInputReader, CommandLineError
from whatshap.vcf import VcfReader, VcfError, VariantTable, VariantCallPhase, VcfInvalidChromosome
from whatshap.core import NumericSampleIds
from whatshap.timer import StageTimer
from whatshap.utils import Region, stdout_is_regular_file


logger = logging.getLogger(__name__)


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    arg("-o", "--output",
        default=sys.stdout,
        help="Output file. If omitted, use standard output.")
    arg("--reference", "-r", metavar="FASTA",
        help="Reference file. Must be accompanied by .fai index (create with samtools faidx)")
    arg("--no-reference", action="store_true", default=False,
        help="Detect alleles without requiring a reference, at the expense of phasing quality "
        "(in particular for long reads)")
    arg("--regions", dest="regions", metavar="REGION", default=None, action="append",
        help="Specify region(s) of interest to limit the tagging to reads/variants "
        "overlapping those regions. You can specify a space-separated list of "
        "regions in the form of chrom:start-end, chrom (consider entire chromosome), "
        "or chrom:start (consider region from this start to end of chromosome).")
    arg("--ignore-linked-read", default=False, action="store_true",
        help="Ignore linkage information stored in BX tags of the reads.")
    arg("--linked-read-distance-cutoff", "-d", metavar="LINKEDREADDISTANCE",
        default=50000, type=int,
        help="Assume reads with identical BX tags belong to different read clouds if their "
        "distance is larger than LINKEDREADDISTANCE (default: %(default)s).")
    arg("--ignore-read-groups", default=False, action="store_true",
        help="Ignore read groups in BAM/CRAM header and assume all reads come "
        "from the same sample.")
    arg("--sample", dest="given_samples", metavar="SAMPLE", default=None, action="append",
        help="Name of a sample to phase. If not given, all samples in the "
        "input VCF are phased. Can be used multiple times.")
    arg("--output-haplotag-list", dest="haplotag_list", metavar="HAPLOTAG_LIST", default=None,
        help="Write assignments of read names to haplotypes (tab separated) to given "
        "output file. If filename ends in .gz, then output is gzipped.")
    arg("--tag-supplementary", default=False, action="store_true",
        help="Also tag supplementary alignments. Supplementary alignments are assigned to the "
        "same haplotype as the primary alignment (default: only tag primary alignments).")
    arg("--ploidy", metavar="PLOIDY", default=2, type=int, help="Ploidy (default: %(default)s).")
    arg("--skip-missing-contigs", default=False, action="store_true",
        help="Skip reads that map to a contig that does not exist in the VCF")
    arg("--output-threads", "--out-threads", default=1, type=int,
        help="Number of threads to use for output file writing (passed to pysam). "
        "For optimal performance, instead pipe output into 'samtools view' to compress.")
    arg("variant_file", metavar="VCF", help="VCF file with phased variants "
        "(must be gzip-compressed and indexed)")
    arg("alignment_file", metavar="ALIGNMENTS",
        help="BAM/CRAM file with alignments to be tagged by haplotype")
# fmt: on


def validate(args, parser):
    if args.reference is not None and args.no_reference:
        parser.error("Options --reference and --no-reference cannot be used together")


def md5_of(filename):
    with open(filename, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def get_variant_information(variant_table: VariantTable, sample: str):
    """
    Return (vpos_to_phase_info, variants) pair, where vpos_to_phase_info maps
    a variant position to a pair (integer block_id, phase0),
    and variants is a list of all non-homozygous variants.
    """
    genotypes = variant_table.genotypes_of(sample)
    phases: List[Optional[VariantCallPhase]] = variant_table.phases_of(sample)

    vpos_to_phase_info = dict()
    variants = []
    for v, gt, phase in zip(variant_table.variants, genotypes, phases):
        if phase is None or phase.block_id is None:
            continue
        # map block_id to tuple of phases
        phase_info = int(phase.block_id), phase.phase
        vpos_to_phase_info[v.position] = phase_info
        if not gt.is_homozygous():
            variants.append(v)

    return vpos_to_phase_info, variants


def attempt_add_phase_information(
    alignment, read_to_haplotype, bxtag_to_haplotype, linked_read_cutoff, ignore_linked_read
):
    is_tagged = 0
    haplotype_name = "none"
    phaseset = "none"
    try:
        haplotype, quality, phaseset = read_to_haplotype[alignment.query_name]
        haplotype_name = f"H{haplotype + 1}"
        alignment.set_tag("HP", haplotype + 1)
        alignment.set_tag("PC", quality)
        alignment.set_tag("PS", phaseset)
        is_tagged = 1
    except KeyError:
        # check if reads with same tag have been assigned
        if not ignore_linked_read:
            try:
                tag = alignment.get_tag("BX")
            except KeyError:
                read_clouds = []
            else:  # alignment has BX tag
                read_clouds = bxtag_to_haplotype[tag]

            for reference_start, haplotype, phaseset in read_clouds:
                if abs(reference_start - alignment.reference_start) <= linked_read_cutoff:
                    haplotype_name = f"H{haplotype + 1}"
                    alignment.set_tag("HP", haplotype + 1)
                    alignment.set_tag("PC", value=None)
                    alignment.set_tag("PS", phaseset)
                    is_tagged = 1
                    break
    return is_tagged, haplotype_name, phaseset


def load_chromosome_variants(
    vcf_reader: VcfReader, chromosome: str, regions: Sequence[Tuple[int, Optional[int]]]
) -> VariantTable:
    try:
        logger.debug(f"Loading variants from {len(regions)} distinct region(s)")
        variant_table = vcf_reader.fetch_regions(chromosome, regions)
        logger.debug(f"Loaded {len(variant_table)} variants for chromosome {chromosome} in VCF")
    except OSError as err:
        # not entirely clear to me why this could raise
        # an OSError at this point?
        logger.error(str(err))
        raise err
    return variant_table


def prepare_haplotag_information(
    variant_table,
    shared_samples,
    phased_input_reader,
    regions,
    ignore_linked_read,
    linked_read_cutoff,
    ploidy,
):
    """
    Read all reads for this chromosome once to create one core.ReadSet per sample.
    This allows to assign phase to paired-end reads based on both reads
    """
    n_multiple_phase_sets = 0
    BX_tag_to_haplotype = defaultdict(list)
    # maps read name to (haplotype, quality, phaseset)
    read_to_haplotype = {}

    for sample in shared_samples:
        variantpos_to_phaseinfo, variants = get_variant_information(variant_table, sample)
        read_set, _ = phased_input_reader.read(
            variant_table.chromosome, variants, sample, regions=regions
        )

        # map BX tag to list of reads
        bx_tag_to_readlist = defaultdict(list)
        if not ignore_linked_read:
            for read in read_set:
                if read.has_BX_tag():
                    bx_tag_to_readlist[read.BX_tag].append(read)

        # all reads processed so far
        processed_reads = set()
        for read in read_set:
            if read.name in processed_reads:
                continue
            # mapping: phaseset --> costs of assigning reads to haplotypes
            haplotype_costs = defaultdict(lambda: [0] * ploidy)

            processed_reads.add(read.name)
            reads_to_consider = {read}

            # reads with same BX tag need to be considered too (unless --ignore-linked-read is set)
            if not ignore_linked_read and read.has_BX_tag():
                for r in bx_tag_to_readlist[read.BX_tag]:
                    if r.name not in processed_reads:
                        # only select reads close to current one
                        if abs(read.reference_start - r.reference_start) <= linked_read_cutoff:
                            reads_to_consider.add(r)

            for r in reads_to_consider:
                processed_reads.add(r.name)
                for v in r:
                    assert v.allele in [0, 1]
                    phaseset, phasing = variantpos_to_phaseinfo[v.position]
                    for hap_index, hap_allele in enumerate(phasing):
                        if v.allele == hap_allele:
                            haplotype_costs[phaseset][hap_index] += v.quality

            l = list(haplotype_costs.items())
            # sort by maximum quality score
            l.sort(key=lambda t: max(t[1]), reverse=True)
            # logger.info('Read %s: %s', read.name, str(l))

            if len(l) == 0:
                continue
            if len(l) > 1:
                n_multiple_phase_sets += 1
            phaseset, scores = l[0]

            # find best and second best haplotype scores for this phaseset
            scores_list = list(enumerate(scores))
            scores_list.sort(key=lambda t: t[1], reverse=True)
            first_ht, first_score = scores_list[0]
            second_ht, second_score = scores_list[1]
            quality = first_score - second_score

            if quality == 0:
                continue

            if not ignore_linked_read and read.has_BX_tag():
                BX_tag_to_haplotype[read.BX_tag].append((read.reference_start, first_ht, phaseset))

            for r in reads_to_consider:
                read_to_haplotype[r.name] = (first_ht, quality, phaseset)
                logger.debug(
                    "Assigned read {} to haplotype {} with a "
                    "quality of {} based on {} covered variants".format(
                        r.name, first_ht, quality, len(r)
                    )
                )
    return BX_tag_to_haplotype, read_to_haplotype, n_multiple_phase_sets


def normalize_user_regions(
    user_regions: Optional[Sequence[str]], bam_references: List[str]
) -> Dict[str, List[Tuple[int, Optional[int]]]]:
    """
    Process and accept user input of the following forms:

    chr:start-end -> chr, start, end
    chr -> chr, 0, None
    chr:start -> chr, start, None

    User input is interpreted as 1-based closed intervals.
    In pysam, coordinates are treated as 0-based half-open,
    so convert user input chr1:1-100 into chr1:0-100

    Args:
        user_regions: list of user input regions
        bam_references: references of BAM file

    Returns:
        dict of lists containing normalized regions per chromosome
    """
    regions: Dict[str, List[Tuple[int, Optional[int]]]] = defaultdict(list)
    if user_regions is None:
        for reference in bam_references:
            regions[reference].append((0, None))
    else:
        bam_references = set(bam_references)
        for region_spec in user_regions:
            region = Region.parse(region_spec)
            if region.chromosome not in bam_references:
                raise ValueError(
                    "Requested reference '{region.chromosome}' not found in input BAM/CRAM"
                )
            regions[region.chromosome].append((region.start, region.end))
    return regions


def compute_variant_file_samples_to_use(vcf_samples, user_given_samples, ignore_read_groups):
    """
    Open variant file and load sample information - check if samples in VCF are compatible
    with user specified list of samples.

    return iterable of VCF samples to use
    """
    samples_in_vcf = set(vcf_samples)
    if len(samples_in_vcf) < 1:
        raise VcfError("No samples detected in VCF file; cannot perform haplotagging")
    logger.info(f"Found {len(samples_in_vcf)} sample(s) in input VCF")
    logger.debug(
        "Found the following samples in input VCF: {}".format(" - ".join(sorted(samples_in_vcf)))
    )

    if ignore_read_groups and user_given_samples is None and len(samples_in_vcf) > 1:
        raise ValueError(
            'When setting "--ignore-read-groups" on '
            "a multi-sample VCF, samples to be used must "
            'be specified via the "--sample" parameter.'
        )

    if user_given_samples is None:
        samples_to_use = samples_in_vcf
    else:
        given_samples = user_given_samples
        missing_samples = set(given_samples) - samples_in_vcf
        if len(missing_samples) > 0:
            raise VcfError(
                "The following samples were specified via the "
                '"--sample" parameter, but are not part of the '
                "input VCF: {}".format(sorted(missing_samples))
            )

        samples_to_use = samples_in_vcf.intersection(given_samples)
        logger.info(f"Keeping {len(samples_to_use)} sample(s) for haplo-tagging")
        logger.debug(
            "Keeping the following samples for haplo-tagging: {}".format(
                " - ".join(sorted(samples_to_use))
            )
        )
    return samples_to_use


def compute_shared_samples(bam_reader, ignore_read_groups, vcf_samples):
    """
    Return final samples to use for haplo-tagging
    """
    read_groups = bam_reader.header.get("RG", [])
    bam_samples = {(rg["SM"] if "SM" in rg else "") for rg in read_groups}

    logger.info(f"Found {len(bam_samples)} sample(s) in BAM file")
    logger.debug(
        "Found the following samples in BAM file: {}".format(",".join(sorted(bam_samples)))
    )

    if not ignore_read_groups:
        shared_samples = bam_samples.intersection(vcf_samples)
        if len(shared_samples) == 0:
            raise ValueError(
                "No common samples between VCF and BAM file detected. "
                'You may restart the analysis setting "--ignore-read-groups" '
                "(if appropriate) to avoid this error."
            )
        elif len(shared_samples) < len(bam_samples):
            missing_samples = " | ".join(sorted(bam_samples - shared_samples))
            logger.warning(
                "Ignoring the following sample(s) for haplo-tagging "
                "because they are not part of the VCF or "
                'were not requested via "--sample": {}'.format(missing_samples)
            )
        else:
            # situation is ok
            pass
    else:
        shared_samples = vcf_samples
    return shared_samples


def open_output_alignment_file(aln_output, reference, vcf_md5, bam_header, threads=1):
    """
    :param aln_output:
    :param reference:
    :param vcf_md5:
    :param bam_header:
    :param exit_stack:
    :return:
    """
    # Prepare header
    # TODO: convince pysam to allow @HS header line
    command_line = " ".join(["whatshap"] + sys.argv[1:])
    PG_entry = {
        "ID": "whatshap",
        "PN": "whatshap",
        "VN": __version__,
        "CL": command_line,
        "m5": vcf_md5,
    }
    if "PG" in bam_header:
        bam_header["PG"].append(PG_entry)
    else:
        bam_header["PG"] = [PG_entry]
    if aln_output is None:
        aln_output = "-"
        kwargs = dict()
    elif str(aln_output).endswith(".cram"):  # FIXME hard-coded value
        if reference is None:
            raise ValueError(
                'Writing CRAM output requires FASTA reference file given via "--reference"'
            )
        kwargs = dict(mode="wc", reference_filename=reference)
    else:
        # Write BAM, disable compression when piping
        if aln_output is sys.stdout and not stdout_is_regular_file():
            kwargs = dict(mode="wb0", threads=threads)
        else:
            kwargs = dict(mode="wb", threads=threads)
    try:
        bam_writer = pysam.AlignmentFile(
            aln_output, header=pysam.AlignmentHeader.from_dict(bam_header), **kwargs
        )
    except OSError as err:
        raise CommandLineError(
            f"Error while initializing alignment output file at path: {aln_output}\n{err}"
        )

    return bam_writer


def open_haplotag_writer(path: str) -> TextIO:
    try:
        writer = xopen(path, mode="wt")
    except OSError as err:
        raise CommandLineError(
            f"Error while initializing haplotag list output at path: {path}\n{err}"
        )
    logger.debug("Writing header line to haplotag list output file")
    print("#readname", "haplotype", "phaseset", "chromosome", sep="\t", file=writer)
    return writer


def ignore_read(alignment, tag_supplementary):
    """
    If supplementary alignments should also be tagged,
    this should only take the haplo-tag of the primary
    alignment into account - this leads to:

    We ignore an alignment [aln]:
    - IF aln is_unmapped OR is_secondary
    - IF tag_supplementary AND aln is_secondary
    - IF not tag_supplementary AND is_supplementary

    :param alignment:
    :param tag_supplementary:
    :return:
    """
    # TODO: could be that some checks here are not needed
    # due to default filtering in ReadSetReader::_usableAlignments

    if alignment.is_unmapped or alignment.is_secondary:
        # unmapped or secondary alignments are never tagged
        ignore = True
    elif tag_supplementary and alignment.is_supplementary:
        # from the previous if, we know
        # the alignment to be primary
        ignore = False
    elif alignment.is_supplementary:
        # tag_supplementary is False, so discard
        ignore = True
    else:
        # whatever is left should be good
        ignore = False
    return ignore


def contigs_with_alignments(af: pysam.AlignmentFile) -> FrozenSet[str]:
    has_alignments = []
    for contig in af.references:
        for _ in af.fetch(contig=contig):
            has_alignments.append(contig)
            break
    return frozenset(has_alignments)


def run_haplotag(
    variant_file,
    alignment_file,
    output=None,
    reference: Union[None, bool, str] = False,
    regions=None,
    ignore_linked_read=False,
    given_samples=None,
    linked_read_distance_cutoff=50000,
    ignore_read_groups: bool = False,
    haplotag_list: Optional[List[str]] = None,
    tag_supplementary: bool = False,
    skip_missing_contigs: bool = False,
    output_threads: int = 1,
    ploidy: int = 2,
):
    timers = StageTimer()
    timers.start("haplotag-run")

    if output in (None, sys.stdout) and sys.stdout.isatty():
        raise CommandLineError(
            "Refusing to write BAM to the terminal. Either use the '-o' option or redirect "
            "standard output with '>'."
        )
    with ExitStack() as stack:
        timers.start("haplotag-init")
        try:
            vcf_reader = stack.enter_context(
                VcfReader(variant_file, only_snvs=False, phases=True, ploidy=ploidy)
            )
        except OSError as err:
            raise CommandLineError(f"Error while loading variant file {variant_file}: {err}")

        use_vcf_samples = compute_variant_file_samples_to_use(
            vcf_reader.samples, given_samples, ignore_read_groups
        )
        try:
            bam_reader = stack.enter_context(
                pysam.AlignmentFile(
                    alignment_file,
                    reference_filename=reference if reference else None,
                    require_index=True,
                )
            )
        except OSError as err:
            raise CommandLineError(f"Error while loading alignment file {alignment_file}: {err}")
        # This checks also sample compatibility with VCF
        shared_samples = compute_shared_samples(bam_reader, ignore_read_groups, use_vcf_samples)

        # Check if user has specified a subset of regions per chromosome
        user_regions = normalize_user_regions(regions, bam_reader.references)

        include_unmapped = regions is None
        phased_input_reader = stack.enter_context(
            PhasedInputReader(
                [alignment_file],
                None if reference is False else reference,
                NumericSampleIds(),
                ignore_read_groups,
                only_snvs=False,
                duplicates=True,
            )
        )
        if phased_input_reader.has_alignments and reference is None:
            raise CommandLineError(
                "A reference FASTA needs to be provided with -r/--reference; "
                "or use --no-reference at the expense of phasing quality."
            )

        bam_writer = stack.enter_context(
            open_output_alignment_file(
                output,
                reference,
                md5_of(variant_file),
                bam_reader.header.to_dict(),
                threads=output_threads,
            )
        )
        if haplotag_list is not None:
            haplotag_writer = stack.enter_context(open_haplotag_writer(haplotag_list))
        else:
            haplotag_writer = None

        timers.stop("haplotag-init")
        logger.debug(
            "All input/output files initialized (time: {})".format(timers.elapsed("haplotag-init"))
        )
        timers.start("haplotag-process")

        n_alignments = 0
        n_tagged = 0
        n_multiple_phase_sets = 0

        has_alignments = contigs_with_alignments(bam_reader)

        for chrom, regions in user_regions.items():
            logger.debug(f"Processing chromosome {chrom}")

            if chrom not in has_alignments:
                # Skip chromosomes without alignments. This allows to have extra chromosomes in the
                # BAM header compared to the VCF.
                continue
            try:
                variant_table = load_chromosome_variants(vcf_reader, chrom, regions)
            except VcfInvalidChromosome:
                if skip_missing_contigs:
                    logger.info(
                        f"Skipping reads on '{chrom}' because the contig does not exist in the VCF"
                    )
                    continue
                else:
                    raise CommandLineError(
                        f"Input BAM/CRAM contains reads on contig '{chrom}', but that contig does "
                        "not exist in the VCF header. To bypass this check, use "
                        "--skip-missing-contigs"
                    )
            except VcfError as e:
                raise CommandLineError(str(e))
            if variant_table is not None:
                logger.debug("Preparing haplotype information")
                (BX_tag_to_haplotype, read_to_haplotype, n_mult) = prepare_haplotag_information(
                    variant_table,
                    shared_samples,
                    phased_input_reader,
                    regions,
                    ignore_linked_read,
                    linked_read_distance_cutoff,
                    ploidy,
                )
                n_multiple_phase_sets += n_mult
            else:
                # avoid uninitialized variables
                BX_tag_to_haplotype = None
                read_to_haplotype = None

            assert not include_unmapped or len(regions) == 1
            for start, end in regions:
                logger.debug("Working on %s:%s-%s", chrom, start, end)
                for alignment in bam_reader.fetch(contig=chrom, start=start, stop=end):
                    n_alignments += 1
                    haplotype_name = "none"
                    phaseset = "none"

                    if variant_table is None or ignore_read(alignment, tag_supplementary):
                        # - If no variants in VCF for this chromosome,
                        # alignments just get written to output
                        # - Ignored reads are simply
                        # written to the output BAM
                        # Existing tags HP, PC and PS are removed
                        alignment.set_tag("HP", value=None)
                        alignment.set_tag("PC", value=None)
                        alignment.set_tag("PS", value=None)
                    else:
                        (is_tagged, haplotype_name, phaseset) = attempt_add_phase_information(
                            alignment,
                            read_to_haplotype,
                            BX_tag_to_haplotype,
                            linked_read_distance_cutoff,
                            ignore_linked_read,
                        )
                        n_tagged += is_tagged

                        if not is_tagged:
                            # Remove any existing tags HP, PC and PS if the aligment does
                            # not have phasing information
                            alignment.set_tag("HP", value=None)
                            alignment.set_tag("PC", value=None)
                            alignment.set_tag("PS", value=None)

                    bam_writer.write(alignment)
                    if haplotag_writer is not None and not (
                        alignment.is_secondary or alignment.is_supplementary
                    ):
                        print(
                            alignment.query_name,
                            haplotype_name,
                            phaseset,
                            chrom,
                            sep="\t",
                            file=haplotag_writer,
                        )

                    if n_alignments % 100_000 == 0:
                        logger.debug(f"Processed {n_alignments} alignment records.")
        if include_unmapped:
            logger.debug("Copying unmapped reads to output")
            for alignment in bam_reader.fetch(until_eof=True):
                bam_writer.write(alignment)
        timers.stop("haplotag-process")
        logger.debug("Processing complete (time: {})".format(timers.elapsed("haplotag-process")))

    timers.stop("haplotag-run")

    logger.info("\n== SUMMARY ==")
    logger.info("Total alignments processed:              %12d", n_alignments)
    logger.info("Alignments that could be tagged:         %12d", n_tagged)
    logger.info("Alignments spanning multiple phase sets: %12d", n_multiple_phase_sets)
    logger.info("Finished in %.1f s", timers.elapsed("haplotag-run"))


def main(args):
    if args.no_reference:
        args.reference = False
    del args.no_reference
    run_haplotag(**vars(args))
