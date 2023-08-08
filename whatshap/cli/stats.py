"""
Print phasing statistics of a single VCF file
"""
import logging
from collections import defaultdict
from contextlib import ExitStack
import dataclasses
from statistics import median
from typing import List, Tuple, Optional, Dict, Sequence, Iterator
from math import isnan

from ..vcf import VcfReader, VcfVariant, VariantTable

logger = logging.getLogger(__name__)


# fmt: off
def add_arguments(parser):
    add = parser.add_argument
    add("--gtf", default=None, help="Write phased blocks to GTF file.")
    add("--sample", metavar="SAMPLE", help="Name of the sample "
        "to process. If not given, use first sample found in VCF.")
    add("--chr-lengths", metavar="FILE",
        help="Override chromosome lengths in VCF with those from FILE (one line per chromosome, "
        "tab separated '<chr> <length>'). Lengths are used to compute NG50 values.")
    add("--tsv", metavar="FILE", help="Write statistics in tab-separated value format to FILE")
    add("--only-snvs", default=False, action="store_true", help="Only process SNVs "
        "and ignore all other variants.")
    add("--block-list", metavar="FILE", help="Write list of all blocks to FILE (one block per line)")
    add("--chromosome", dest="chromosomes", metavar="CHROMOSOME", default=[], action="append",
        help="Name of chromosome(s) to process. If not given, all chromosomes in the "
        "input VCF are considered. Can be used multiple times and accepts a comma-separated list. ")
    add("vcf", metavar="VCF", help="Phased VCF file")
# fmt: on


def validate(args, parser):
    pass


class PhasedBlock:
    def __init__(self, chromosome=None):
        self.phases = {}
        self.leftmost_variant = None
        self.rightmost_variant = None
        self.chromosome = chromosome

    def add(self, variant, phase):
        if len(self.phases) == 0:
            self.leftmost_variant = variant
            self.rightmost_variant = variant
        else:
            if variant < self.leftmost_variant:
                self.leftmost_variant = variant
            if self.rightmost_variant < variant:
                self.rightmost_variant = variant
        self.phases[variant] = phase

    def span(self):
        """Returns the length of the covered genomic region in bp."""
        return self.rightmost_variant.position - self.leftmost_variant.position

    def variants(self):
        return list(sorted(self.phases.keys()))

    def count_snvs(self):
        return sum(int(variant.is_snv()) for variant in self.phases)

    def split(self, split_left: int, split_right: int) -> Tuple["PhasedBlock", "PhasedBlock"]:
        """Split this phaseblock in two, based on given positions. The first phaseblock will contain
        the variants to the left of split_left and the second the variants to the right of split_right."""
        assert split_left <= split_right
        left_block = PhasedBlock(chromosome=self.chromosome)
        right_block = PhasedBlock(chromosome=self.chromosome)
        for variant, phase in self.phases.items():
            if variant.position < split_left:
                left_block.add(variant, phase)
            elif variant.position > split_right:
                right_block.add(variant, phase)
        return left_block, right_block

    def __repr__(self):
        return f"PhasedBlock({str(self.phases)})"

    def __len__(self):
        return len(self.phases)

    def __lt__(self, other):
        return (self.leftmost_variant, self.rightmost_variant) < (
            other.leftmost_variant,
            other.rightmost_variant,
        )


class GtfWriter:
    def __init__(self, file):
        self._file = file

    def write(self, chromosome, start, stop, name):
        """
        Write a feature to the GTF. start is 0-based.
        """
        assert start < stop
        print(
            chromosome,
            "Phasing",
            "exon",
            start + 1,
            stop,
            ".",
            "+",
            ".",
            f'gene_id "{name}"; transcript_id "{name}.1";',
            sep="\t",
            file=self._file,
        )


@dataclasses.dataclass
class DetailedStats:
    variants: int = 0
    phased: int = 0
    unphased: int = 0
    singletons: int = 0
    blocks: int = 0
    variant_per_block_median: float = float("nan")
    variant_per_block_avg: float = float("nan")
    variant_per_block_min: int = 0
    variant_per_block_max: int = 0
    variant_per_block_sum: int = 0
    bp_per_block_median: float = float("nan")
    bp_per_block_avg: float = float("nan")
    bp_per_block_min: int = 0
    bp_per_block_max: int = 0
    bp_per_block_sum: int = 0
    heterozygous_variants: int = 0
    heterozygous_snvs: int = 0
    phased_snvs: int = 0
    phased_fraction: float = 0.0
    phased_snvs_fraction: float = 0.0
    block_n50: float = float("nan")

    def print(self):
        # Parameters for value formatting
        max_integer_width = max(
            len(str(int(value))) for value in vars(self).values() if not isnan(value)
        )
        value_width = max(max_integer_width, 8)
        format_int = f"{value_width}.0f"
        format_float = f"{value_width + 3}.2f"
        format_percent = f"{value_width + 3}.1%"
        format_param = ">21"

        # fmt: off
        print(
            f"{'Variants in VCF':{format_param}}: {self.variants:{format_int}}",
            f"{'Heterozygous':{format_param}}: {self.heterozygous_variants:{format_int}}    ({self.heterozygous_snvs:{format_int}}    SNVs)",
            f"{'Phased':{format_param}}: {self.phased:{format_int}}    ({self.phased_snvs:{format_int}}    SNVs)",
            f"{'Heterozygous phased':{format_param}}: {self.phased_fraction:{format_percent}} ({self.phased_snvs_fraction:{format_percent}} SNVs)",
            f"{'Unphased':{format_param}}: {self.unphased:{format_int}}    (not considered below)",
            f"{'Singletons':{format_param}}: {self.singletons:{format_int}}    (not considered below)",
            f"{'Blocks':{format_param}}: {self.blocks:{format_int}}",
            "",
            "Block sizes (no. of variants)",
            f"{'Sum of sizes':{format_param}}: {self.variant_per_block_sum:{format_int}}    variants",
            f"{'Median block size':{format_param}}: {self.variant_per_block_median:{format_float}} variants",
            f"{'Average block size':{format_param}}: {self.variant_per_block_avg:{format_float}} variants",
            f"{'Largest block':{format_param}}: {self.variant_per_block_max:{format_int}}    variants",
            f"{'Smallest block':{format_param}}: {self.variant_per_block_min:{format_int}}    variants",
            "",
            "Block lengths (basepairs)",
            f"{'Sum of lengths':{format_param}}: {self.bp_per_block_sum:{format_int}}    bp",
            f"{'Median block length':{format_param}}: {self.bp_per_block_median:{format_float}} bp",
            f"{'Average block length':{format_param}}: {self.bp_per_block_avg:{format_float}} bp",
            f"{'Longest block':{format_param}}: {self.bp_per_block_max:{format_int}}    bp",
            f"{'Shortest block':{format_param}}: {self.bp_per_block_min:{format_int}}    bp",
            f"{'Block NG50':{format_param}}: {self.block_n50:{format_int}}    bp",
            sep="\n"
        )
        # fmt: on
        assert self.phased + self.unphased + self.singletons == self.heterozygous_variants


def n50(lengths: List[int], target_length: Optional[int] = None) -> int:
    if target_length is None:
        target_length = sum(lengths)

    lengths.sort(reverse=True)
    total = 0
    for length in lengths:
        total += length
        if total >= 0.5 * target_length:
            return length
    return 0


def compute_ng50(blocks: List[PhasedBlock], chr_lengths: Dict[str, int]):
    chromosomes = {b.chromosome for b in blocks}
    target_length = 0
    for chromosome in sorted(chromosomes):
        try:
            target_length += chr_lengths[chromosome]
        except KeyError:
            logger.warning(
                "Not able to compute NG50 because length of contig '%s' not available", chromosome
            )
            return float("nan")

    block_lengths = [b.span() for b in blocks]
    return n50(block_lengths, target_length=target_length)


class PhasingStats:
    def __init__(self):
        self.blocks = []
        self.split_blocks = []
        self.unphased = 0
        self.variants = 0
        self.heterozygous_variants = 0
        self.heterozygous_snvs = 0
        self.phased_snvs = 0

    def __iadd__(self, other):
        self.blocks.extend(other.blocks)
        self.split_blocks.extend(other.split_blocks)
        self.unphased += other.unphased
        self.variants += other.variants
        self.heterozygous_variants += other.heterozygous_variants
        self.heterozygous_snvs += other.heterozygous_snvs
        self.phased_snvs += other.phased_snvs
        return self

    def add_blocks(self, blocks: Sequence[PhasedBlock]):
        self.blocks.extend(blocks)
        self.split_blocks.extend(self.get_nonoverlapping_blocks())

    def add_unphased(self, unphased: int = 1):
        self.unphased += unphased

    def add_variants(self, variants: int):
        self.variants += variants

    def add_heterozygous_variants(self, variants: int):
        self.heterozygous_variants += variants

    def add_heterozygous_snvs(self, snvs: int):
        self.heterozygous_snvs += snvs

    def get_nonoverlapping_blocks(self) -> List[PhasedBlock]:
        """Split phase blocks into nonoverlapping subblocks"""
        pos_sorted_blocks = sorted(
            self.blocks, key=lambda b: (b.chromosome, b.leftmost_variant.position), reverse=True
        )

        # filter out blocks with only one variant
        pos_sorted_blocks = [b for b in pos_sorted_blocks if len(b) > 1]

        # iterate over blocks and split if overlapping until no blocks remain.
        split_blocks = []
        while pos_sorted_blocks:
            block = pos_sorted_blocks.pop()
            if pos_sorted_blocks:
                block_end = block.rightmost_variant.position
                next_block = pos_sorted_blocks[-1]
                next_block_start = next_block.leftmost_variant.position
                next_block_end = next_block.rightmost_variant.position

                # Check if next block overlapps current. If so split the current block.
                if (block_end > next_block_start) and (block.chromosome == next_block.chromosome):
                    block, new_block = block.split(next_block_start, next_block_end)

                    # Update sorting if right-side block is added.
                    if len(new_block) > 1:
                        pos_sorted_blocks.append(new_block)
                        pos_sorted_blocks = sorted(
                            pos_sorted_blocks,
                            key=lambda b: (b.chromosome, b.leftmost_variant.position),
                            reverse=True,
                        )

                    # Skip the left-side block if is it too short after splitting
                    if len(block) < 2:
                        continue
            split_blocks.append(block)

        return split_blocks

    def get_detailed_stats(self, chr_lengths: Optional[Dict[str, int]] = None) -> DetailedStats:
        """Return DetailedStats"""
        block_sizes = sorted(len(block) for block in self.blocks if len(block) > 1)
        n_singletons = sum(1 for block in self.blocks if len(block) == 1)
        # Block length stats calculated from split interleaved blocks to avoid inflating values
        block_lengths = sorted(block.span() for block in self.split_blocks if len(block) > 1)
        phased_snvs = sum(block.count_snvs() for block in self.blocks if len(block) > 1)
        if block_sizes:
            return DetailedStats(
                variants=self.variants,
                phased=sum(block_sizes),
                unphased=self.unphased,
                singletons=n_singletons,
                blocks=len(block_sizes),
                variant_per_block_median=median(block_sizes),
                variant_per_block_avg=sum(block_sizes) / len(block_sizes)
                if len(block_sizes)
                else float("nan"),
                variant_per_block_min=block_sizes[0],
                variant_per_block_max=block_sizes[-1],
                variant_per_block_sum=sum(block_sizes),
                bp_per_block_median=median(block_lengths),
                bp_per_block_avg=sum(block_lengths) / len(block_lengths)
                if len(block_lengths)
                else float("nan"),
                bp_per_block_min=block_lengths[0],
                bp_per_block_max=block_lengths[-1],
                bp_per_block_sum=sum(block_lengths),
                heterozygous_variants=self.heterozygous_variants,
                heterozygous_snvs=self.heterozygous_snvs,
                phased_snvs=phased_snvs,
                phased_fraction=sum(block_sizes) / self.heterozygous_variants
                if self.heterozygous_variants
                else float("nan"),
                phased_snvs_fraction=phased_snvs / self.heterozygous_snvs
                if self.heterozygous_snvs
                else float("nan"),
                block_n50=compute_ng50(self.split_blocks, chr_lengths)
                if chr_lengths is not None
                else float("nan"),
            )
        else:
            return DetailedStats(
                variants=self.variants,
                unphased=self.unphased,
                singletons=n_singletons,
                heterozygous_variants=self.heterozygous_variants,
                heterozygous_snvs=self.heterozygous_snvs,
            )


def unpack_chromosomes(chromosomes: List[str]) -> List[str]:
    """Unpack list chromosomes by splitting comma-separated entries."""
    unpacked = (chromosome for entry in chromosomes for chromosome in entry.split(","))
    return [chromosome for chromosome in unpacked if chromosome != ""]


def parse_chr_lengths(filename) -> Dict[str, int]:
    """
    Parse chromosome lengths from file filename. The file should have two columns with
    chromosome names and lengths respectively.
    """
    chr_lengths = {}
    with open(filename) as f:
        for line in f:
            fields = line.split("\t")
            assert len(fields) == 2
            chr_lengths[fields[0]] = int(fields[1])
    return chr_lengths


def parse_variant_tables(
    vcf_reader: VcfReader, chromosomes: Optional[Sequence[str]] = None
) -> Iterator[VariantTable]:
    """
    Parse variant_tables from vcf_reader. If chromosomes are given and VCF is indexed,
    theses are accessed by direct lookup.
    """
    if chromosomes and vcf_reader.index_exists():
        for chromosome in chromosomes:
            yield vcf_reader.fetch(chromosome)
    else:
        yield from vcf_reader


def get_chr_lengths(
    vcf_reader: VcfReader, chr_lengths_file: Optional[str] = None
) -> Dict[str, int]:
    """
    Return a dictionary that maps a chromosome name to the chromosome’s length. The
    mapping is read from chr_lengths_file if provided, and from the VCF header otherwise.
    """
    if chr_lengths_file:
        chr_lengths = parse_chr_lengths(chr_lengths_file)
        logger.info("Read length of %d chromosomes from %s", len(chr_lengths), chr_lengths_file)
    else:
        chr_lengths = {
            contig.name: contig.length
            for contig in vcf_reader.contigs.values()
            if contig.length is not None
        }
        if not chr_lengths:
            logger.warning(
                "VCF header does not contain contig lengths, cannot compute NG50. "
                "Consider using --chr-lengths"
            )
    return chr_lengths


def write_to_block_list(
    block_list_file, blocks: Dict[int, PhasedBlock], chromosome: str, sample: str
):
    """
    Write phase blocks for chromosome to block_list_file.
    """
    block_ids = sorted(blocks.keys())
    for block_id in block_ids:
        print(
            sample,
            chromosome,
            block_id,
            blocks[block_id].leftmost_variant.position + 1,
            blocks[block_id].rightmost_variant.position + 1,
            len(blocks[block_id]),
            sep="\t",
            file=block_list_file,
        )


@dataclasses.dataclass
class GtfBlock:
    start: Optional[int] = 0
    end: Optional[int] = 0
    id: Optional[int] = None

    def add(self, variant: VcfVariant):
        self.end = variant.position + 1


def get_phase_blocks(
    chromosome: str,
    gtfwriter: GtfWriter,
    sample: str,
    stats: PhasingStats,
    variant_table: VariantTable,
) -> Dict[int, PhasedBlock]:
    """
    Parse phase blocks from variant_table for sample. Returns map of block ids to phaseblocks.
    """
    genotypes = variant_table.genotypes_of(sample)
    phases = variant_table.phases_of(sample)
    assert len(genotypes) == len(phases) == len(variant_table.variants)

    blocks: Dict[int, PhasedBlock] = defaultdict(PhasedBlock)
    prev_block = GtfBlock()
    for variant, genotype, phase in zip(variant_table.variants, genotypes, phases):
        stats.add_variants(1)
        if genotype.is_homozygous():
            continue
        stats.add_heterozygous_variants(1)
        if variant.is_snv():
            stats.add_heterozygous_snvs(1)

        if phase is None:
            stats.add_unphased()
            continue

        blocks[phase.block_id].add(variant, phase)
        if gtfwriter:
            if prev_block.id is None:
                prev_block = GtfBlock(variant.position, variant.position + 1, phase.block_id)
            else:
                if prev_block.id != phase.block_id:
                    gtfwriter.write(chromosome, prev_block.start, prev_block.end, prev_block.id)
                    prev_block = GtfBlock(variant.position, variant.position + 1, phase.block_id)

                prev_block.add(variant)

    # Add chromosome information to each block. This is needed to
    # sort blocks later when we compute NG50s
    for block_id, block in blocks.items():
        block.chromosome = chromosome

    if gtfwriter and prev_block.id is not None:
        gtfwriter.write(chromosome, prev_block.start, prev_block.end, prev_block.id)

    return blocks


def run_stats(
    vcf,
    sample=None,
    gtf=None,
    tsv=None,
    block_list=None,
    only_snvs=False,
    chromosomes=None,
    chr_lengths=None,
):
    gtfwriter = tsv_file = block_list_file = None

    if chromosomes is not None:
        chromosomes = unpack_chromosomes(chromosomes)

    with ExitStack() as stack:
        if gtf:
            gtf_file = stack.enter_context(open(gtf, "wt"))
            gtfwriter = GtfWriter(gtf_file)

        vcf_reader = VcfReader(vcf, phases=True, only_snvs=only_snvs)
        if len(vcf_reader.samples) == 0:
            logger.error("Input VCF does not contain any sample")
            return 1
        else:
            logger.info(f"Found {len(vcf_reader.samples)} sample(s) in input VCF")
        if sample:
            if sample in vcf_reader.samples:
                sample = sample
            else:
                logger.error(f"Requested sample ({sample}) not found")
                return 1
        else:
            sample = vcf_reader.samples[0]
            logger.info(f"Reporting results for sample {sample}")

        chr_lengths = get_chr_lengths(vcf_reader, chr_lengths)

        if tsv:
            tsv_file = stack.enter_context(open(tsv, "w"))
            field_names = [f.name for f in dataclasses.fields(DetailedStats)]
            print("#sample", "chromosome", "file_name", *field_names, sep="\t", file=tsv_file)

        if block_list:
            block_list_file = stack.enter_context(open(block_list, "w"))
            print(
                "#sample",
                "chromosome",
                "phase_set",
                "from",
                "to",
                "variants",
                sep="\t",
                file=block_list_file,
            )

        print(f"Phasing statistics for sample {sample} from file {vcf}")
        total_stats = PhasingStats()
        given_chromosomes = chromosomes
        seen_chromosomes = set()
        for variant_table in parse_variant_tables(vcf_reader, given_chromosomes):
            chromosome = variant_table.chromosome
            seen_chromosomes.add(chromosome)
            if given_chromosomes and chromosome not in given_chromosomes:
                continue

            stats = PhasingStats()
            print(f"---------------- Chromosome {chromosome} ----------------")
            blocks = get_phase_blocks(chromosome, gtfwriter, sample, stats, variant_table)

            if block_list_file:
                write_to_block_list(block_list_file, blocks, chromosome, sample)

            stats.add_blocks(blocks.values())

            detailed_stats = stats.get_detailed_stats(chr_lengths)
            detailed_stats.print()
            if tsv_file:
                print(sample, chromosome, vcf, sep="\t", end="\t", file=tsv_file)
                print(*dataclasses.astuple(detailed_stats), sep="\t", file=tsv_file)

            total_stats += stats

            if given_chromosomes and set(given_chromosomes) <= seen_chromosomes:
                break

        if len(seen_chromosomes) > 1:
            print("---------------- ALL chromosomes (aggregated) ----------------")
            detailed_stats = total_stats.get_detailed_stats(chr_lengths)
            detailed_stats.print()
            if tsv_file:
                print(sample, "ALL", vcf, sep="\t", end="\t", file=tsv_file)
                print(*dataclasses.astuple(detailed_stats), sep="\t", file=tsv_file)


def main(args):
    run_stats(**vars(args))
