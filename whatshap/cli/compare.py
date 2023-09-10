"""
Compare two or more phased variant files
"""
import logging
import math
from collections import defaultdict
from contextlib import ExitStack
import dataclasses
from itertools import chain, permutations
from typing import Set, List, Optional, DefaultDict, Dict

from whatshap.vcf import VcfReader, VcfVariant, VariantTable, PloidyError
from whatshap.core import Genotype
from whatshap.polyphase.solver import SwitchFlipCalculator
from whatshap.cli import CommandLineError


logger = logging.getLogger(__name__)

COUNT_WIDTH = 9


# fmt: off
def add_arguments(parser):
    add = parser.add_argument
    add('--sample', metavar='SAMPLE', default=None, help='Name of the sample '
        'to process. If not given, use first sample found in VCF.')
    add('--names', metavar='NAMES', default=None, help='Comma-separated list '
        'of data set names to be used in the report (in same order as VCFs).')
    add('--ignore-sample-name', default=False, action='store_true', help='For single '
        'sample VCFs, ignore sample name and assume all samples are the same.')
    add('--tsv-pairwise', metavar='TSVPAIRWISE', default=None, help='Filename to write '
        'comparison results from pair-wise comparison to (tab-separated).')
    add('--tsv-multiway', metavar='TSVMULTIWAY', default=None, help='Filename to write '
        'comparison results from multiway comparison to (tab-separated). Only for diploid VCFs.')
    add('--only-snvs', default=False, action="store_true", help='Only process SNVs '
        'and ignore all other variants.')
    add('--switch-error-bed', default=None, help='Write BED file with switch error positions '
        'to given filename. Only for diploid VCFs.')
    add('--plot-blocksizes', default=None, help='Write PDF file with a block length histogram '
        'to given filename (requires matplotlib).')
    add('--plot-sum-of-blocksizes', default=None, help='Write PDF file with a block length histogram in which the height of each bar corresponds to the sum of lengths.')
    add('--longest-block-tsv', default=None, help='Write position-wise agreement of longest '
        'joint blocks in each chromosome to tab-separated file. Only for diploid VCFs.')
    add('--ploidy', '-p', metavar='PLOIDY', type=int, default=2, help='The ploidy of the sample(s) (default: %(default)s).')
    # TODO: what's the best way to request "two or more" VCFs?
    add('vcf', nargs='+', metavar='VCF/BCF', help='At least two phased variant files (VCF or BCF) to be compared.')
# fmt: on


def validate(args, parser):
    if len(args.vcf) < 2:
        parser.error("At least two VCFs need to be given.")
    if args.ploidy < 2:
        parser.error("Ploidy must be > 1.")
    if args.ploidy > 2 and args.tsv_multiway:
        parser.error("Option --tsv-multiway can only be used if ploidy=2.")
    if args.ploidy > 2 and args.switch_error_bed:
        parser.error("Option --switch-error-bed can only be used if ploidy=2.")
    if args.ploidy > 2 and args.longest_block_tsv:
        parser.error("Option --longest-block-tsv can only be used if ploidy=2.")


class SwitchFlips:
    def __init__(self, switches: int = 0, flips: int = 0):
        self.switches: int = switches
        self.flips: int = flips

    def __iadd__(self, other):
        self.switches += other.switches
        self.flips += other.flips
        return self

    def __repr__(self):
        return f"SwitchFlips(switches={self.switches}, flips={self.flips})"

    def __str__(self):
        return f"{self.switches}/{self.flips}"


class PhasingErrors:
    def __init__(
        self,
        switches: int = 0,
        hamming: int = 0,
        switch_flips: Optional[SwitchFlips] = None,
        diff_genotypes: int = 0,
    ):
        self.switches = switches
        self.hamming = hamming
        self.switch_flips = SwitchFlips() if switch_flips is None else switch_flips
        self.diff_genotypes = diff_genotypes

    def __iadd__(self, other: object) -> "PhasingErrors":
        if not isinstance(other, PhasingErrors):
            raise TypeError("Can only add to PhasingErrors")
        self.switches += other.switches
        self.hamming += other.hamming
        self.switch_flips += other.switch_flips
        self.diff_genotypes += other.diff_genotypes
        return self

    def __repr__(self):
        return "PhasingErrors(switches={}, hamming={}, switch_flips={}, diff_genotypes={})".format(
            self.switches, self.hamming, self.switch_flips, self.diff_genotypes
        )


def complement(s):
    """
    >>> complement('01100')
    '10011'
    """
    t = {"0": "1", "1": "0"}
    return "".join(t[c] for c in s)


def hamming(s0, s1):
    """
    >>> hamming('ABCD', 'AXCY')
    2
    """
    assert len(s0) == len(s1)
    return sum(c0 != c1 for c0, c1 in zip(s0, s1))


def switch_encoding(phasing):
    """
    >>> switch_encoding('0001011')
    '001110'
    """
    assert isinstance(phasing, str)
    return "".join(("0" if phasing[i - 1] == phasing[i] else "1") for i in range(1, len(phasing)))


def compute_switch_flips(phasing0, phasing1) -> SwitchFlips:
    """
    >>> compute_switch_flips("00011", "00100")
    SwitchFlips(switches=1, flips=0)
    >>> compute_switch_flips("00011", "00111")
    SwitchFlips(switches=0, flips=1)
    >>> compute_switch_flips("000", "001")
    SwitchFlips(switches=1, flips=0)
    """
    assert len(phasing0) == len(phasing1)
    s0 = switch_encoding(phasing0)
    s1 = switch_encoding(phasing1)
    result = SwitchFlips()
    switches_in_a_row = 0
    for i, (p0, p1) in enumerate(zip(s0, s1)):
        if p0 != p1:
            switches_in_a_row += 1
        if (i + 1 == len(s0)) or (p0 == p1):
            result.flips += switches_in_a_row // 2
            result.switches += switches_in_a_row % 2
            switches_in_a_row = 0

    return result


def compute_matching_genotype_pos(phasing0, phasing1):
    """
    Computes the positions on which both phasings agree on the genotype.
    """
    assert len(phasing0) == len(phasing1)
    assert len(phasing0) >= 2
    assert len(phasing0[0]) == len(phasing1[0])
    assert all(len(phasing0[i]) == len(phasing0[0]) for i in range(1, len(phasing0)))
    num_vars = len(phasing0[0])
    matching_pos = [
        i
        for i in range(num_vars)
        if Genotype([int(hap[i]) for hap in phasing0])
        == Genotype([int(hap[i]) for hap in phasing1])
    ]
    return matching_pos


def compute_switch_errors_poly(phasing0, phasing1, matching_pos=None):
    """
    Computes the number of necessary switches to transform phasing 0 into phasing 1 or vice versa.
    Positions with non-matching genotypes are omitted.
    """
    assert len(phasing0) == len(phasing1)
    assert len(phasing0) >= 2
    assert len(phasing0[0]) == len(phasing1[0])
    assert all(len(phasing0[i]) == len(phasing0[0]) for i in range(1, len(phasing0)))
    num_vars = len(phasing0[0])

    # If positions with matching genotypes are not precomputed, do it here!
    if matching_pos is None:
        matching_pos = compute_matching_genotype_pos(phasing0, phasing1)

    phasing0_matched = ["".join([hap[i] for i in matching_pos]) for hap in phasing0]
    phasing1_matched = ["".join([hap[i] for i in matching_pos]) for hap in phasing1]

    vector_error = compute_switch_flips_poly(
        phasing0_matched,
        phasing1_matched,
        switch_cost=1,
        flip_cost=2 * num_vars * len(phasing0) + 1,
    )
    assert vector_error.flips == 0

    return vector_error.switches


def compute_switch_flips_poly(phasing0, phasing1, switch_cost=1, flip_cost=1):
    """
    Computes the combined number of switches and flips, which are needed to transform phasing 0 into
    phasing 1 or vice versa.
    """
    (result, switches_in_column, flips_in_column, poswise_config) = compute_switch_flips_poly_bt(
        phasing0, phasing1, switch_cost=switch_cost, flip_cost=flip_cost
    )
    return result


def compute_switch_flips_poly_bt(
    phasing0, phasing1, report_error_positions=False, switch_cost=1, flip_cost=1
):
    # Check input
    if len(phasing0) != len(phasing1):
        logger.error(
            "Incompatible phasings. Number of haplotypes is not equal "
            f"({len(phasing0)} != {len(phasing1)})."
        )
    assert len(phasing0) == len(phasing1)

    num_pos = len(phasing0[0])
    if num_pos == 0:
        return SwitchFlips(), None, None, None
    ploidy = len(phasing0)
    if ploidy == 0:
        return SwitchFlips(), None, None, None
    for i in range(0, len(phasing1)):
        if len(phasing1[i]) != num_pos:
            logger.error(
                "Inconsistent input for phasing. Haplotypes have different lengths "
                f"( len(phasing1[0]={num_pos} != len(phasing1[{i}]={len(phasing1[i])}."
            )
        assert len(phasing1[i]) == num_pos
        if len(phasing0[i]) != num_pos:
            logger.error(
                "Inconsistent input for phasing. Haplotypes have different lengths "
                f"( len(phasing1[0]={num_pos} != len(phasing0[{i}]={len(phasing0[i])}."
            )
        assert len(phasing1[i]) == num_pos
    if ploidy > 6:
        logger.warning(
            "Computing vector error with more than 6 haplotypes. This may take very long ..."
        )

    # Compute comparison
    calc = SwitchFlipCalculator(ploidy, switch_cost, flip_cost)
    result = SwitchFlips()
    (
        switches,
        flips,
        switches_in_column,
        flips_in_column,
        positionwise_config,
    ) = calc.compute_switch_flips_poly(phasing0, phasing1)

    # Aggregate results
    result.switches = switches / ploidy
    result.flips = flips / ploidy
    return result, switches_in_column, flips_in_column, positionwise_config


def poly_num_switches(perm0, perm1):
    cost = 0
    for i in range(len(perm0)):
        if perm0[i] != perm1[i]:
            cost += 1
    return cost


def compare_block(phasing0, phasing1):
    """Input are two lists of haplotype sequences over {0,1}."""
    assert len(phasing0) == len(phasing1)
    ploidy = len(phasing0)

    minimum_hamming_distance = float("inf")
    # compute minimum hamming distance
    for permutation in permutations(phasing0):
        # compute sum of hamming distances
        total_hamming = 0
        for i in range(ploidy):
            total_hamming += hamming(phasing1[i], permutation[i])
        total_hamming /= float(ploidy)
        minimum_hamming_distance = min(minimum_hamming_distance, total_hamming)

    matching_pos = compute_matching_genotype_pos(phasing0, phasing1)

    if ploidy == 2:
        # conversion to int is allowed, as there should be no fractional error counts for diploid comparisons
        switches = int(hamming(switch_encoding(phasing0[0]), switch_encoding(phasing1[0])))
        switch_flips = compute_switch_flips(phasing0[0], phasing1[0])
        minimum_hamming_distance = int(minimum_hamming_distance)
    else:
        switches = compute_switch_errors_poly(phasing0, phasing1, matching_pos)
        switch_flips = compute_switch_flips_poly(phasing0, phasing1)

    return PhasingErrors(
        switches=switches,
        hamming=minimum_hamming_distance,
        switch_flips=switch_flips,
        diff_genotypes=len(phasing0[0]) - len(matching_pos),
    )


def fraction2percentstr(nominator, denominator):
    if denominator == 0:
        return "--"
    else:
        return f"{nominator * 100.0 / denominator:.2f}%"


def safefraction(nominator, denominator):
    if denominator == 0:
        return float("nan")
    else:
        return nominator / denominator


class BedCreator:
    def __init__(self, chromosome: str, dataset_names: List[str]):
        self._chromosome = chromosome
        self._annotation = "{}<-->{}".format(*dataset_names)

    def records(self, phasing0, phasing1, positions):
        """
        Determine positions of switch errors between two phasings
        and yield one BED record per switch error (encoded as a tuple).
        The annotation_string is added to each record.
        """
        assert len(phasing0) == len(phasing1) == len(positions)
        switch_encoding0 = switch_encoding(phasing0)
        switch_encoding1 = switch_encoding(phasing1)
        for i, (sw0, sw1) in enumerate(zip(switch_encoding0, switch_encoding1)):
            if sw0 != sw1:
                yield (self._chromosome, positions[i] + 1, positions[i + 1] + 1, self._annotation)


def print_stat(text: str, value=None, value2=None, text_width=37):
    """
    Print a line like this:

         text: value
    """
    text = text.rjust(text_width)
    if value is None:
        assert value2 is None
        print(text)
    else:
        if value == "-":
            value = "-" * COUNT_WIDTH
        else:
            value = str(value).rjust(COUNT_WIDTH)
        if value2 is None:
            print(text + ":", value)
        else:
            print(text + ":", value, str(value2).rjust(COUNT_WIDTH))


def print_errors(errors, phased_pairs):
    print_stat("phased pairs of variants assessed", phased_pairs)
    print_stat("switch errors", errors.switches)
    print_stat("switch error rate", fraction2percentstr(errors.switches, phased_pairs))
    print_stat("switch/flip decomposition", errors.switch_flips)
    print_stat(
        "switch/flip rate",
        fraction2percentstr(errors.switch_flips.switches + errors.switch_flips.flips, phased_pairs),
    )


@dataclasses.dataclass
class PairwiseComparisonResults:
    intersection_blocks: int
    covered_variants: int
    all_assessed_pairs: int
    all_switches: int
    all_switch_rate: float
    all_switchflips: SwitchFlips
    all_switchflip_rate: float
    blockwise_hamming: int
    blockwise_hamming_rate: int
    blockwise_diff_genotypes: int
    blockwise_diff_genotypes_rate: int
    largestblock_assessed_pairs: int
    largestblock_switches: int
    largestblock_switch_rate: float
    largestblock_switchflips: SwitchFlips
    largestblock_switchflip_rate: float
    largestblock_hamming: int
    largestblock_hamming_rate: float
    largestblock_diff_genotypes: int
    largestblock_diff_genotypes_rate: float


@dataclasses.dataclass
class BlockStats:
    variant_count: int
    span: int


def collect_common_variants(
    variant_tables: List[VariantTable], sample_names: List[str]
) -> Set[VcfVariant]:
    common_variants = None
    for variant_table, sample in zip(variant_tables, sample_names):
        het_variants = [
            v
            for v, gt in zip(variant_table.variants, variant_table.genotypes_of(sample))
            if not gt.is_homozygous()
        ]
        if common_variants is None:
            common_variants = set(het_variants)
        else:
            common_variants.intersection_update(het_variants)
    assert common_variants is not None
    return common_variants


def compare(
    variant_tables: List[VariantTable],
    sample_names: List[str],
    dataset_names: List[str],
    ploidy: int,
):
    """
    Return a PairwiseComparisonResults object if the variant_tables has a length of 2.
    """
    assert len(variant_tables) > 1

    common_variants = collect_common_variants(variant_tables, sample_names)
    assert common_variants is not None

    print_stat("common heterozygous variants", len(common_variants))
    print_stat("(restricting to these below)")
    phases = []
    sorted_variants = sorted(common_variants, key=lambda v: v.position)
    for variant_table, sample in zip(variant_tables, sample_names):
        p = [
            phase
            for variant, phase in zip(variant_table.variants, variant_table.phases_of(sample))
            if variant in common_variants
        ]
        assert [v for v in variant_table.variants if v in common_variants] == sorted_variants
        assert len(p) == len(common_variants)
        phases.append(p)

    # blocks[variant_table_index][block_id] is a list of indices into common_variants
    blocks: List[DefaultDict[int, List[int]]] = [defaultdict(list) for _ in variant_tables]
    block_intersection = defaultdict(list)
    for variant_index in range(len(common_variants)):
        any_none = False
        for i in range(len(phases)):
            phase = phases[i][variant_index]
            if phase is None or any(p is None for p in phase.phase):
                any_none = True
            else:
                blocks[i][phase.block_id].append(variant_index)
        if not any_none:
            joint_block_id = tuple(
                phase[variant_index].block_id for phase in phases  # type: ignore
            )
            block_intersection[joint_block_id].append(variant_index)

    # create statistics on each block in each data set
    block_stats = compute_block_stats(blocks, sorted_variants)

    for dataset_name, blck in zip(dataset_names, blocks):
        print_stat(
            f"non-singleton blocks in {dataset_name}",
            len([b for b in blck.values() if len(b) > 1]),
        )
        print_stat("--> covered variants", sum(len(b) for b in blck.values() if len(b) > 1))

    intersection_block_count = sum(1 for b in block_intersection.values() if len(b) > 1)
    intersection_block_variants = sum(len(b) for b in block_intersection.values() if len(b) > 1)
    print_stat("non-singleton intersection blocks", intersection_block_count)
    print_stat("--> covered variants", intersection_block_variants)
    if len(variant_tables) == 2:
        (
            bed_records,
            longest_block_agreement,
            longest_block_positions,
            pairwise_comparison,
        ) = compare_pair(
            block_intersection,
            intersection_block_count,
            intersection_block_variants,
            phases,
            ploidy,
            sorted_variants,
            BedCreator(variant_tables[0].chromosome, dataset_names),
        )

        return (
            pairwise_comparison,
            bed_records,
            block_stats,
            longest_block_positions,
            longest_block_agreement,
            None,
        )
    else:
        assert ploidy == 2
        multiway_results = compare_multiway(block_intersection, dataset_names, phases)
        return None, None, block_stats, None, None, multiway_results


def compare_pair(
    block_intersection,
    intersection_block_count,
    intersection_block_variants,
    phases,
    ploidy,
    sorted_variants,
    bed_creator: Optional[BedCreator],
):
    longest_block = 0
    longest_block_errors = PhasingErrors()
    longest_block_positions = []
    longest_block_agreement = []
    phased_pairs = 0
    bed_records = []
    total_errors = PhasingErrors()
    total_compared_variants = 0
    for block in block_intersection.values():
        if len(block) < 2:
            continue
        phasing0 = []
        phasing1 = []
        for j in range(ploidy):
            p0 = "".join(str(phases[0][i].phase[j]) for i in block)
            p1 = "".join(str(phases[1][i].phase[j]) for i in block)
            phasing0.append(p0)
            phasing1.append(p1)
        block_positions = [sorted_variants[i].position for i in block]
        errors = compare_block(phasing0, phasing1)

        # TODO: extend to polyploid
        if ploidy == 2 and bed_creator is not None:
            bed_records.extend(bed_creator.records(phasing0[0], phasing1[0], block_positions))
        total_errors += errors
        phased_pairs += len(block) - 1
        total_compared_variants += len(block)
        if len(block) > longest_block:
            longest_block = len(block)
            longest_block_errors = errors
            longest_block_positions = block_positions
            # TODO: extend to polyploid
            if ploidy == 2:
                if hamming(phasing0, phasing1) < hamming(phasing0[0], complement(phasing1[0])):
                    longest_block_agreement = [
                        1 * (p0 == p1) for p0, p1 in zip(phasing0[0], phasing1[0])
                    ]
                else:
                    longest_block_agreement = [
                        1 * (p0 != p1) for p0, p1 in zip(phasing0[0], phasing1[0])
                    ]
    longest_block_assessed_pairs = max(longest_block - 1, 0)
    print_stat("ALL INTERSECTION BLOCKS", "-")
    print_errors(total_errors, phased_pairs)
    print_stat("Block-wise Hamming distance", total_errors.hamming)
    print_stat(
        "Block-wise Hamming distance [%]",
        fraction2percentstr(total_errors.hamming, total_compared_variants),
    )
    print_stat("Different genotypes", total_errors.diff_genotypes)
    print_stat(
        "Different genotypes [%]",
        fraction2percentstr(total_errors.diff_genotypes, total_compared_variants),
    )
    print_stat("LARGEST INTERSECTION BLOCK", "-")
    print_errors(longest_block_errors, longest_block_assessed_pairs)
    print_stat("Hamming distance", longest_block_errors.hamming)
    print_stat(
        "Hamming distance [%]", fraction2percentstr(longest_block_errors.hamming, longest_block)
    )
    print_stat("Different genotypes", longest_block_errors.diff_genotypes)
    print_stat(
        "Different genotypes [%]",
        fraction2percentstr(longest_block_errors.diff_genotypes, longest_block),
    )
    pcr = PairwiseComparisonResults(
        intersection_blocks=intersection_block_count,
        covered_variants=intersection_block_variants,
        all_assessed_pairs=phased_pairs,
        all_switches=total_errors.switches,
        all_switch_rate=safefraction(total_errors.switches, phased_pairs),
        all_switchflips=total_errors.switch_flips,
        all_switchflip_rate=safefraction(
            total_errors.switch_flips.switches + total_errors.switch_flips.flips, phased_pairs
        ),
        blockwise_hamming=total_errors.hamming,
        blockwise_hamming_rate=safefraction(total_errors.hamming, total_compared_variants),
        blockwise_diff_genotypes=total_errors.diff_genotypes,
        blockwise_diff_genotypes_rate=safefraction(
            total_errors.diff_genotypes, total_compared_variants
        ),
        largestblock_assessed_pairs=longest_block_assessed_pairs,
        largestblock_switches=longest_block_errors.switches,
        largestblock_switch_rate=safefraction(
            longest_block_errors.switches, longest_block_assessed_pairs
        ),
        largestblock_switchflips=longest_block_errors.switch_flips,
        largestblock_switchflip_rate=safefraction(
            longest_block_errors.switch_flips.switches + longest_block_errors.switch_flips.flips,
            longest_block_assessed_pairs,
        ),
        largestblock_hamming=longest_block_errors.hamming,
        largestblock_hamming_rate=safefraction(longest_block_errors.hamming, longest_block),
        largestblock_diff_genotypes=longest_block_errors.diff_genotypes,
        largestblock_diff_genotypes_rate=safefraction(
            longest_block_errors.diff_genotypes, longest_block
        ),
    )
    return bed_records, longest_block_agreement, longest_block_positions, pcr


def compare_multiway(block_intersection, dataset_names, phases):
    histogram = defaultdict(int)
    total_compared = 0
    for block in block_intersection.values():
        if len(block) < 2:
            continue
        total_compared += len(block) - 1
        phasings = ["".join(str(phases[j][i].phase[0]) for i in block) for j in range(len(phases))]
        switch_encodings = [switch_encoding(p) for p in phasings]
        for i in range(len(block) - 1):
            s = "".join(switch_encodings[j][i] for j in range(len(switch_encodings)))
            s = min(s, complement(s))
            histogram[s] += 1
    print_stat("Compared pairs of variants", total_compared)
    bipartitions = list(histogram.keys())
    bipartitions.sort()
    multiway_results = {}  # (dataset_list0, dataset_list1) --> count
    for i, s in enumerate(bipartitions):
        count = histogram[s]
        if i == 0:
            assert {c for c in s} == set("0")
            print("ALL AGREE")
        elif i == 1:
            print("DISAGREEMENT")
        left, right = [], []
        for name, leftright in zip(dataset_names, s):
            if leftright == "0":
                left.append(name)
            else:
                right.append(name)
        print_stat(
            ("{{{}}} vs. {{{}}}".format(",".join(left), ",".join(right))),
            count,
            fraction2percentstr(count, total_compared),
        )
        multiway_results[(",".join(left), ",".join(right))] = count
    return multiway_results


def compute_block_stats(
    blocks: List[DefaultDict[int, List[int]]], sorted_variants: List[VcfVariant]
):
    block_stats = []
    for block in blocks:
        l = []
        for block_id, variant_indices in block.items():
            if len(variant_indices) < 2:
                continue
            span = (
                sorted_variants[variant_indices[-1]].position
                - sorted_variants[variant_indices[0]].position
            )
            l.append(BlockStats(len(variant_indices), span))
        block_stats.append(l)
    return block_stats


def create_blocksize_histogram(filename, block_stats, names, use_weights=False):
    try:
        import matplotlib
        import numpy

        matplotlib.use("pdf")
        from matplotlib import pyplot
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError:
        raise CommandLineError(
            "To use option --plot-blocksizes, you need to have numpy and matplotlib installed."
        )

    assert len(block_stats) == len(names)

    color_list = ["#ffa347", "#0064c8", "#b42222", "#22a5b4", "#b47c22", "#6db6ff"]
    if len(color_list) < len(block_stats):
        color_count = len(block_stats)
        color_list = pyplot.cm.Set1([n / color_count for n in range(color_count)])
    colors = color_list[: len(block_stats)]

    with PdfPages(filename) as pdf:
        for what, xlabel in [
            (lambda stats: stats.variant_count, "variant count"),
            (lambda stats: stats.span, "span [bp]"),
        ]:
            pyplot.figure(figsize=(10, 8))
            max_value = max(what(stats) for stats in chain(*block_stats))
            common_bins = numpy.logspace(0, math.ceil(math.log10(max_value)), 50)
            for l, name, color in zip(block_stats, names, colors):
                x = [what(stats) for stats in l]
                n, bins, patches = pyplot.hist(
                    x,
                    bins=common_bins,
                    alpha=0.6,
                    color=color,
                    label=name,
                    weights=x if use_weights else None,
                )
            pyplot.xlabel(xlabel)
            pyplot.ylabel("Number of blocks")
            pyplot.gca().set_xscale("log")
            pyplot.gca().set_yscale("log")
            pyplot.grid(True)
            pyplot.legend()
            pdf.savefig()
            pyplot.close()

            pyplot.figure(figsize=(10, 8))
            common_bins = numpy.logspace(0, math.ceil(math.log10(max_value)), 25)
            x = [[what(stats) for stats in l] for l in block_stats]
            n, bins, patches = pyplot.hist(
                x,
                bins=common_bins,
                alpha=0.6,
                color=colors,
                label=names,
                weights=x if use_weights else None,
            )
            pyplot.xlabel(xlabel)
            pyplot.ylabel("Number of blocks")
            pyplot.gca().set_xscale("log")
            pyplot.gca().set_yscale("log")
            pyplot.grid(True)
            pyplot.legend()
            pdf.savefig()
            pyplot.close()


def run_compare(
    vcf,
    ploidy,
    names=None,
    sample=None,
    ignore_sample_name=False,
    tsv_pairwise=None,
    tsv_multiway=None,
    only_snvs=False,
    switch_error_bed=None,
    plot_blocksizes=None,
    plot_sum_of_blocksizes=None,
    longest_block_tsv=None,
):
    vcf_readers = [
        VcfReader(f, only_snvs=only_snvs, phases=True, ploidy=ploidy, mav=True) for f in vcf
    ]
    if names:
        dataset_names = names.split(",")
        if len(dataset_names) != len(vcf):
            raise CommandLineError(
                "Number of names given with --names does not equal number of VCFs."
            )
    else:
        dataset_names = [f"file{i}" for i in range(len(vcf))]

    sample_names = get_sample_names(
        vcf_readers, requested_sample=sample, ignore_name=ignore_sample_name
    )

    with ExitStack() as stack:
        tsv_pairwise_file = tsv_multiway_file = longest_block_tsv_file = switch_error_bedfile = None
        if tsv_pairwise:
            tsv_pairwise_file = stack.enter_context(open(tsv_pairwise, "w"))

        if tsv_multiway:
            tsv_multiway_file = stack.enter_context(open(tsv_multiway, "w"))
            print(
                "#sample",
                "chromosome",
                "dataset_list0",
                "dataset_list1",
                "count",
                sep="\t",
                file=tsv_multiway_file,
            )

        if longest_block_tsv:
            longest_block_tsv_file = stack.enter_context(open(longest_block_tsv, "w"))
            print(
                "#dataset_name0",
                "dataset_name1",
                "#sample",
                "chromosome",
                "position",
                "phase_agreeing",
                sep="\t",
                file=longest_block_tsv_file,
            )

        if tsv_pairwise_file:
            fields = [
                "#sample",
                "chromosome",
                "dataset_name0",
                "dataset_name1",
                "file_name0",
                "file_name1",
            ]
            field_names = [f.name for f in dataclasses.fields(PairwiseComparisonResults)]
            fields.extend(field_names)
            fields.extend(["het_variants0", "only_snvs"])
            print(*fields, sep="\t", file=tsv_pairwise_file)

        if switch_error_bed:
            switch_error_bedfile = stack.enter_context(open(switch_error_bed, "w"))

        if len(set(sample_names)) > 1 and ignore_sample_name:
            print(
                "Comparing phasings for samples:",
                ", ".join(sample_names),
                " (--ignore-sample-names selected)",
            )
        else:
            print("Comparing phasings for sample", sample_names[0])

        vcfs = get_variant_tables(vcf_readers, vcf)
        chromosomes = get_common_chromosomes(vcfs)
        if len(chromosomes) == 0:
            raise CommandLineError("No chromosome is contained in all VCFs. Aborting.")
        logger.info("Chromosomes present in all VCFs: %s", ", ".join(chromosomes))

        print("FILENAMES")
        longest_name = max(len(n) for n in dataset_names)
        for name, filename in zip(dataset_names, vcf):
            print(name.rjust(longest_name + 2), "=", filename)

        width = max(longest_name, 15) + 5

        all_block_stats = [[] for _ in vcfs]

        def add_block_stats(block_stats):
            assert len(block_stats) == len(all_block_stats)
            for big_list, new_list in zip(all_block_stats, block_stats):
                big_list.extend(new_list)

        for chromosome in sorted(chromosomes):
            print(f"---------------- Chromosome {chromosome} ----------------")
            all_bed_records = []
            variant_tables = [vcf[chromosome] for vcf in vcfs]
            all_variants_union = set()
            all_variants_intersection = None
            het_variants_union = set()
            het_variants_intersection = None
            het_variant_sets = []
            het_variants0 = None
            print("VARIANT COUNTS (heterozygous / all): ")
            for variant_table, name, sample in zip(variant_tables, dataset_names, sample_names):
                all_variants_union.update(variant_table.variants)
                het_variants = [
                    v
                    for v, gt in zip(variant_table.variants, variant_table.genotypes_of(sample))
                    if not gt.is_homozygous()
                ]
                if het_variants0 is None:
                    het_variants0 = len(het_variants)
                het_variants_union.update(het_variants)
                if all_variants_intersection is None:
                    all_variants_intersection = set(variant_table.variants)
                    het_variants_intersection = set(het_variants)
                else:
                    all_variants_intersection.intersection_update(variant_table.variants)
                    het_variants_intersection.intersection_update(het_variants)
                het_variant_sets.append(set(het_variants))
                print(
                    f"{name}:".rjust(width),
                    str(len(het_variants)).rjust(COUNT_WIDTH),
                    "/",
                    str(len(variant_table.variants)).rjust(COUNT_WIDTH),
                )
            print(
                "UNION:".rjust(width),
                str(len(het_variants_union)).rjust(COUNT_WIDTH),
                "/",
                str(len(all_variants_union)).rjust(COUNT_WIDTH),
            )
            print(
                "INTERSECTION:".rjust(width),
                str(len(het_variants_intersection)).rjust(COUNT_WIDTH),
                "/",
                str(len(all_variants_intersection)).rjust(COUNT_WIDTH),
            )

            for i in range(len(vcfs)):
                for j in range(i + 1, len(vcfs)):
                    print(
                        "PAIRWISE COMPARISON: {} <--> {}:".format(
                            dataset_names[i], dataset_names[j]
                        )
                    )
                    (
                        results,
                        bed_records,
                        block_stats,
                        longest_block_positions,
                        longest_block_agreement,
                        multiway_results,
                    ) = compare(
                        [variant_tables[i], variant_tables[j]],
                        [sample_names[i], sample_names[j]],
                        [dataset_names[i], dataset_names[j]],
                        ploidy,
                    )
                    if len(vcfs) == 2:
                        add_block_stats(block_stats)
                    all_bed_records.extend(bed_records)
                    sample_name = (
                        f"{sample_names[i]}_{sample_names[j]}"
                        if ignore_sample_name
                        else sample_names[i]
                    )
                    if tsv_pairwise_file:
                        fields = [
                            sample_name,
                            chromosome,
                            dataset_names[i],
                            dataset_names[j],
                            vcf[i],
                            vcf[j],
                        ]
                        fields.extend(dataclasses.astuple(results))
                        fields.extend([het_variants0, int(only_snvs)])
                        print(*fields, sep="\t", file=tsv_pairwise_file)
                    if longest_block_tsv_file:
                        assert ploidy == 2
                        assert len(longest_block_positions) == len(longest_block_agreement)
                        for position, phase_agreeing in zip(
                            longest_block_positions, longest_block_agreement
                        ):
                            print(
                                dataset_names[i],
                                dataset_names[j],
                                sample_name,
                                chromosome,
                                position,
                                phase_agreeing,
                                sep="\t",
                                file=longest_block_tsv_file,
                            )

            # if requested, write all switch errors found in the current chromosome to the bed file
            if switch_error_bedfile:
                assert ploidy == 2
                all_bed_records.sort()
                for record in all_bed_records:
                    print(*record, sep="\t", file=switch_error_bedfile)

            if len(vcfs) > 2:
                assert ploidy == 2
                print("MULTIWAY COMPARISON OF ALL PHASINGS:")
                (
                    results,
                    bed_records,
                    block_stats,
                    longest_block_positions,
                    longest_block_agreement,
                    multiway_results,
                ) = compare(variant_tables, sample_names, dataset_names, ploidy)
                add_block_stats(block_stats)
                if tsv_multiway_file:
                    sample_name = (
                        "_".join(set(sample_names)) if ignore_sample_name else sample_names[0]
                    )
                    for (dataset_list0, dataset_list1), count in multiway_results.items():
                        print(
                            sample_name,
                            chromosome,
                            "{" + dataset_list0 + "}",
                            "{" + dataset_list1 + "}",
                            count,
                            sep="\t",
                            file=tsv_multiway_file,
                        )

        if plot_blocksizes:
            create_blocksize_histogram(plot_blocksizes, all_block_stats, dataset_names)
        if plot_sum_of_blocksizes:
            create_blocksize_histogram(
                plot_sum_of_blocksizes, all_block_stats, dataset_names, use_weights=True
            )


def get_common_chromosomes(vcfs: List[Dict[str, VariantTable]]) -> List[str]:
    common = None
    for chrom_variant_table_map in vcfs:
        chromosomes = chrom_variant_table_map.keys()
        if common is None:
            common = set(chromosomes)
        else:
            common.intersection_update(chromosomes)
    if common is None:
        return []
    return sorted(common)


def get_variant_tables(
    vcf_readers: List[VcfReader], vcf_filenames: List[str]
) -> List[Dict[str, VariantTable]]:
    vcfs = []
    for reader, filename in zip(vcf_readers, vcf_filenames):
        # create dict mapping chromosome names to VariantTables
        m = dict()
        logger.info("Reading phasing from %r", filename)
        try:
            for variant_table in reader:
                m[variant_table.chromosome] = variant_table
        except PloidyError as e:
            raise CommandLineError(f"Provided ploidy is invalid: {e}. Aborting.")
        vcfs.append(m)
    return vcfs


def get_sample_names(
    vcf_readers: List[VcfReader], requested_sample: Optional[str], ignore_name: bool = False
) -> List[str]:
    first_samples = []
    sample_intersection = None
    for vcf_reader in vcf_readers:
        if sample_intersection is None:
            sample_intersection = set(vcf_reader.samples)
        else:
            sample_intersection.intersection_update(vcf_reader.samples)

        if ignore_name and len(vcf_reader.samples) > 1:
            raise CommandLineError(
                "File '{file}' contains multiple samples, option --ignore-sample-name not available.".format(
                    file=vcf_reader.path
                )
            )
        first_samples.append(vcf_reader.samples[0])
    assert sample_intersection is not None
    if requested_sample:
        sample_intersection.intersection_update([requested_sample])
        if len(sample_intersection) == 0:
            raise CommandLineError(
                "Sample {!r} requested on command-line not found in all VCFs".format(
                    requested_sample
                )
            )
        sample_names = [requested_sample] * len(vcf_readers)
    elif ignore_name:
        sample_names = first_samples
    else:
        if len(sample_intersection) == 0:
            raise CommandLineError("None of the samples is present in all VCFs")
        elif len(sample_intersection) == 1:
            sample_names = [list(sample_intersection)[0]] * len(vcf_readers)
        else:
            raise CommandLineError(
                "More than one sample is present in all VCFs, please use"
                " --sample to specify which sample to work on."
            )
    return sample_names


def main(args):
    run_compare(**vars(args))
