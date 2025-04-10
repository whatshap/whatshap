"""
Phase variants in VCF based on information from haplotagged reads
"""

from collections import defaultdict
import itertools
import logging
import sys
from contextlib import ExitStack
from typing import List, Optional, Sequence, Union, Dict, Tuple


from whatshap import __version__
from whatshap.cli import (
    PhasedInputReader,
    CommandLineError,
    log_memory_usage,
    raise_if_any_sample_not_in_vcf,
)
from whatshap.core import NumericSampleIds, Variant, Read
from whatshap.timer import StageTimer
from whatshap.utils import ChromosomeFilter, IndexedFasta
from whatshap.vcf import VcfReader, PhasedVcfWriter, VcfError, VcfVariant, VariantCallPhase

logger = logging.getLogger(__name__)


# fmt: off
def add_arguments(parser):
    arg = parser.add_argument
    arg("-o", "--output",
        default=sys.stdout,
        help="Output file. If omitted, use standard output.")
    arg("--reference", "-r", metavar="FASTA",
        help="Reference file. Must be accompanied by .fai index (create with samtools faidx)")
    arg("--gap-threshold", "-g", metavar="PERCENT", default=70, type=int,
        help="Threshold percentage for qualities. If the percentage of votes for the variant is less than this value, "
        "the algorithm does not assign any information to the variant.")
    arg("--cut-poly", "-c", metavar="LENGTH", default=10, type=int,
        help="Ignore variants within homopolymers longer than the cut value.")
    arg("--only-indels", "-i", default=False, action="store_true",
        help="Add phasing information only to indels.")
    arg("--sample", dest="samples", metavar="SAMPLE", default=[], action="append",
        help="Name of a sample to phase. If not given, all samples in the "
        "input VCF are phased. Can be used multiple times.")
    arg("--ignore-read-groups", default=False, action="store_true",
        help="Ignore read groups in BAM/CRAM header and assume all reads come from the same sample.")
    arg("--chromosome", dest="chromosomes", metavar="CHROMOSOME", default=[], action="append",
        help="Name of chromosome to phase. If not given, all chromosomes in the input VCF are phased. "
        "Can be used multiple times.")
    arg("--no-mav", dest="mav", default=True, action="store_false", help="Ignore multiallelic variants.")
    arg("--exclude-chromosome", dest="excluded_chromosomes", default=[], action="append",
        help="Name of chromosome not to phase.")
    arg("variant_file", metavar="VCF", help="VCF file with variants to phase (must be gzip-compressed and indexed)")
    arg("alignment_file", metavar="ALIGNMENTS",
        help="BAM/CRAM file with alignments tagged by haplotype and phase set")
# fmt: on


def run_haplotagphase(
    variant_file,
    alignment_file,
    output=None,
    samples: Optional[Sequence[str]] = None,
    reference: Union[None, bool, str] = False,
    ignore_read_groups: bool = False,
    only_indels: bool = False,
    chromosomes: Optional[List[str]] = None,
    excluded_chromosomes: Optional[List[str]] = None,
    gap_threshold: int = 70,
    cut_poly: int = 10,
    write_command_line_header: bool = True,
    mav: bool = True,
    tag: str = "PS",
):
    if samples is None:
        samples = []
    if reference is None:
        raise CommandLineError("Option --reference should be specified")
    timers = StageTimer()
    timers.start("haplotagphase-run")
    command_line: Optional[str]
    if write_command_line_header:
        command_line = "(whatshap {}) {}".format(__version__, " ".join(sys.argv[1:]))
    else:
        command_line = None
    with ExitStack() as stack:
        phased_input_reader = stack.enter_context(
            PhasedInputReader(
                [alignment_file],
                None if reference is False else reference,
                NumericSampleIds(),
                ignore_read_groups,
                only_snvs=False,
            )
        )
        try:
            vcf_writer = stack.enter_context(
                PhasedVcfWriter(
                    command_line=command_line,
                    in_path=variant_file,
                    out_file=output,
                    tag=tag,
                    mav=mav,
                )
            )
        except (OSError, VcfError) as e:
            raise CommandLineError(e)

        vcf_reader = stack.enter_context(VcfReader(variant_file, phases=True, mav=mav))

        if ignore_read_groups and not samples and len(vcf_reader.samples) > 1:
            raise CommandLineError(
                "When using --ignore-read-groups on a VCF with "
                "multiple samples, --sample must also be used."
            )

        if not samples:
            samples = vcf_reader.samples

        assert samples is not None

        raise_if_any_sample_not_in_vcf(vcf_reader, samples)

        with timers("read-fasta"):
            fasta = stack.enter_context(IndexedFasta(reference))
        included_chromosomes = ChromosomeFilter(chromosomes, excluded_chromosomes)
        for variant_table in timers.iterate("parse-vcf", vcf_reader):
            chromosome = variant_table.chromosome
            fasta_chr = fasta[chromosome]
            logger.info(f"Processing chromosome {chromosome}...")
            if chromosome not in included_chromosomes:
                logger.info(f"Leaving chromosome {chromosome} unchanged")
                with timers("write-vcf"):
                    vcf_writer.write_unchanged(chromosome)
                continue
            sample_to_super_reads, sample_to_components = (dict(), dict())
            for sample in vcf_reader.samples:
                logger.info(f"Processing sample {sample}")
                genotypes = variant_table.genotypes_of(sample)
                with timers("read-bam"):
                    reads, _ = phased_input_reader.read(
                        chromosome, variant_table.variants, sample, restricted_genotypes=genotypes
                    )
                phases = variant_table.phases_of(sample)
                if sample not in samples:
                    logger.info(f"Skipping sample {sample}")
                    continue
                homozygous = dict()
                change = dict()
                phased = dict()
                # mapping of detected variants to 0/1 and reversed mappings.
                allele_to_id = defaultdict(dict)
                id_to_allele = defaultdict(dict)
                homozygous_number = 0
                phased_number = 0
                for variant, (phase, genotype) in zip(
                    variant_table.variants, zip(phases, genotypes)
                ):
                    for i, v in enumerate(genotype.as_vector()):
                        allele_to_id[variant.position][v] = i
                        id_to_allele[variant.position][i] = v
                    homozygous[variant.position] = genotype.is_homozygous()
                    phased[variant.position] = phase
                    phased_number += phase is not None
                    homozygous_number += genotype.is_homozygous()
                    change[variant.position] = variant
                logger.info(f"Number of homozygous variants is {homozygous_number}")
                logger.info(f"Number of already phased variants is {phased_number}")
                with timers("compute-votes"):
                    votes = compute_votes(homozygous, reads, allele_to_id)
                with timers("compute-consensus"):
                    sample_to_super_reads[sample], sample_to_components[sample] = consensus(
                        only_indels,
                        gap_threshold,
                        cut_poly,
                        fasta_chr,
                        change,
                        phased,
                        votes,
                        id_to_allele,
                    )
            with timers("write-vcf"):
                vcf_writer.write(chromosome, sample_to_super_reads, sample_to_components)
    timers.stop("haplotagphase-run")
    log_time_and_memory_usage(timers)


def log_time_and_memory_usage(timers):
    logger.info("\n# Resource usage")
    log_memory_usage()
    # fmt: off
    logger.info("Finished in :                              %6.1f s", timers.elapsed("haplotagphase-run"))
    logger.info("Time spent reading reference:              %6.1f s", timers.elapsed("read-fasta"))
    logger.info("Time spent reading VCF:                    %6.1f s", timers.elapsed("parse-vcf"))
    logger.info("Time spent writing VCF:                    %6.1f s", timers.elapsed("write-vcf"))
    logger.info("Time spent reading BAM:                    %6.1f s", timers.elapsed("read-bam"))
    logger.info("Time spent computing votes:                %6.1f s", timers.elapsed("compute-votes"))
    logger.info("Time spent spent computing consensus:      %6.1f s", timers.elapsed("compute-consensus"))
    # fmt: on


def consensus(
    only_indels: bool,
    gap_threshold: int,
    cut_homopolymers: int,
    refseq: str,
    change: Dict[int, VcfVariant],
    phased: Dict[int, Optional[VariantCallPhase]],
    votes: Dict[int, Dict[Tuple[int, int], int]],
    id_to_allele: Dict[int, Dict[int, int]],
) -> Tuple[List[List[Read]], Dict[int, int]]:
    """
    Compute a consensus based on voting and filtering criteria.

    This function processes variant votes to create two consensus sequences (super reads),
    taking into account phasing information, gap threshold, and homopolymer's cutoff
    length.

    Args:
        only_indels: If True, only consider non-SNVs for inclusion in the consensus. SNVs are ignored.
        gap_threshold: The minimum percentage of votes a variant must have to be included.
        cut_homopolymers: The cutoff length for homopolymers. Variants within homopolymers longer than this
            length are excluded. A value of <=0 disables this filter.
        reference_chr: The chromosome sequence from a reference FASTA.
        change: A dictionary mapping variant positions to variant.
        phased: A dictionary indicating the phasing status of variants.
            Variants with `None` are considered unphased.
        votes: A dictionary of variant positions to their votes.
            Each vote includes alleles and their corresponding quality scores.
        id_to_allele: A dictionary mapping variant id and positions to the actual allele.

    Returns:
        A tuple containing two elements:
            - super_reads: Two lists of `Variant` objects representing the haplotypes.
            - components: A dictionary representing the ps.

    """
    super_reads = [[], []]
    components = dict()

    for pos, vote in votes.items():
        best_allele, phase_set, fraction, score = best_candidate(vote)
        components[pos] = phase_set
        if phased[pos] is None:
            if 100 * fraction < gap_threshold:
                continue
            if only_indels and change[pos].is_snv():
                continue
            if cut_homopolymers > 0:
                max_length = max(
                    length_of_homopolymer(refseq, pos + 1, 1, cut_homopolymers),
                    length_of_homopolymer(refseq, pos, -1, cut_homopolymers),
                )
                if max_length > cut_homopolymers:
                    continue
        super_reads[0].append(Variant(pos, allele=id_to_allele[pos][best_allele], quality=score))
        super_reads[1].append(
            Variant(pos, allele=id_to_allele[pos][1 - best_allele], quality=score)
        )
    for read in super_reads:
        read.sort(key=lambda x: x.position)
    return super_reads, components


def best_candidate(var: Dict[Tuple[int, int], int]) -> Tuple[int, int, float, int]:
    """
    Compute the proportion of the best candidate's score relative to the total score of all candidates
    and return this score with a candidate.

    Args:
        var: A dictionary of candidate components, where keys are tuples containing a component identifier
            and an additional identifier, and values are the scores.

    Returns:
        Tuple containing four elements:
            - The allele associated with the best candidate component.
            - The phase set of the candidate
            - The quotient of the best candidate's score divided by the total score of all candidates, representing
              the relative significance of the best candidate's score.
            - The score of the best candidate.

    Examples:
        >>> best_candidate({(1, 2): 50, (2, 3): 100, (3, 4): 75})
        (3, 2, 0.4444444444444444, 100)
        >>> best_candidate({(1, 1): 10, (2, 2): 20, (3, 3): 30, (4, 4): 40})
        (4, 4, 0.4, 40)
        >>> best_candidate({(0, 0): 2})
        (0, 0, 1.0, 2)
        >>> best_candidate({(1, 2): 100, (2, 2): 100, (3, 3): 100})
        (2, 1, 0.3333333333333333, 100)
        >>> best_candidate({(5, 5): 200, (6, 6): 300, (7, 7): 500})
        (7, 7, 0.5, 500)
        >>> best_candidate({(1, 2): 50, (2, 3): 100, (3, 4): 75})
        (3, 2, 0.4444444444444444, 100)
        >>> best_candidate({(1, 1): 10, (2, 2): 20, (3, 3): 30, (4, 4): 40})
        (4, 4, 0.4, 40)
    """
    lst = list(var.items())
    lst.sort(key=lambda x: x[-1], reverse=True)
    (phase_set, allele), score = lst[0]
    total = sum(e[-1] for e in lst)
    q = score / total
    return allele, phase_set, q, score


def length_of_homopolymer(ref: str, start: int, step: int, threshold: int) -> int:
    """
    Compute the length of a homopolymer in a reference string.

    Args:
        ref: The reference string.
        start: The starting index in `ref` for the homopolymer sequence.
        step: The step size to use when moving through `ref`.
            This can be used to control the direction and step length for counting
            (e.g., a step of 1 for forward, -1 for backward).
        threshold: The maximum length to count up to. If the count of
            consecutive repeating characters reaches this threshold,
            the counting stops.

    Returns:
        The length of the polymer, which is the count of consecutive repeating
        characters from the start position, not exceeding the threshold.
    Examples:
        >>> length_of_homopolymer("AAABBBCCC", 0, 1, 10)
        3
        >>> length_of_homopolymer("AAABBBCCC", 2, -1, 10)
        3
        >>> length_of_homopolymer("AAABBBCCC", 3, 1, 2)
        2
        >>> length_of_homopolymer("A", 0, 1, 1)
        1
        >>> length_of_homopolymer("AABBBCCCC", 5, 1, 5)
        4
        >>> length_of_homopolymer("", 0, 1, 10)
        0
    """
    res = 0
    for i in itertools.count(start, step):
        if res < threshold and 0 <= i < len(ref) and ref[i] == ref[start]:
            res += 1
        else:
            break
    return res


def compute_votes(
    is_homozygous: Dict[int, bool], reads: List[Read], allele_to_id: Dict[int, Dict[int, int]]
) -> Dict[int, Dict[Tuple[int, int], int]]:
    """
    Compute votes for variants based on read information.

    This function processes a list of reads to accumulate votes for heterozygous
    variant positions. Each read contributes to the vote of a variant based on its
    phasing set (PS) tag, haplotype (HP) tag, and the allele quality. Reads with
    invalid HP tags or those representing more than a diploid configuration are
    skipped and logged. The vote for each variant is a count of its observed alleles
    weighted by their quality, differentiated by the read's phasing information.

    Parameters:
        is_homozygous: A dictionary indicating whether a variant position is homozygous.
        reads: A list of Read objects, each containing information about variants
            observed in the read, including PS and HP tags, variant position, allele, and quality.
        allele_to_id: A dictionary mapping allele positions and an id of an allele to 0/1 indices.

    Returns:
        A dictionary where keys are variant positions and
        values are dictionaries. Each inner dictionary maps a tuple of (phasing set index, haplotype) to
        the total quality score accumulated for that variant.
    """
    votes = dict()
    number_of_skipped = 0
    for read in reads:
        ps, ht = read.PS_tag - 1, read.HP_tag - 1
        if ht < 0 or ps < 0:
            continue
        if ht > 1:
            number_of_skipped += 1
            continue
        for variant in read:
            if is_homozygous[variant.position]:
                continue
            if variant.position not in votes:
                votes[variant.position] = dict()
            if (ps, 0) not in votes[variant.position]:
                votes[variant.position][(ps, 0)] = 0
                votes[variant.position][(ps, 1)] = 0
            votes[variant.position][
                (ps, ht ^ allele_to_id[variant.position][variant.allele])
            ] += variant.quality
    if number_of_skipped > 0:
        logger.warning(
            f"{number_of_skipped} reads were skipped due incorrect HP. The haplotagphase command supports only a diploid input"
        )
    return votes


def main(args):
    run_haplotagphase(**vars(args))
