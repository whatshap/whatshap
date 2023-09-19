"""
Functions for reading VCFs.
"""
import os
import sys
import math
import logging
import itertools
from dataclasses import dataclass
from abc import ABC, abstractmethod
from os import PathLike
from typing import List, Sequence, Dict, Set, Tuple, Iterable, Optional, Union, TextIO, Iterator

from pysam import VariantFile, VariantHeader, VariantRecord
from pysam.libcbcf import VariantRecordSample

from .core import (
    Read,
    ReadSet,
    PhredGenotypeLikelihoods,
    Genotype,
    binomial_coefficient,
    get_max_genotype_ploidy,
    get_max_genotype_alleles,
)
from .utils import warn_once

logger = logging.getLogger(__name__)


class VcfError(Exception):
    pass


class VcfNotSortedError(VcfError):
    pass


class PloidyError(VcfError):
    pass


class VcfIndexMissing(VcfError):
    pass


class VcfInvalidChromosome(VcfError):
    pass


class VcfInvalidAllele(VcfError):
    pass


@dataclass
class VariantCallPhase:
    block_id: int  # numeric id of the phased block
    phase: Tuple[Optional[int], ...]  # alleles representing the phasing. (1, 0) is 1|0
    quality: Optional[int]


class VcfVariant(ABC):
    """A variant in a VCF file (not to be confused with core.Variant)"""

    position: int
    reference_allele: str

    @abstractmethod
    def get_ref_allele(self):
        pass

    @abstractmethod
    def get_alt_allele(self):
        pass

    @abstractmethod
    def get_alt_allele_list(self):
        pass

    @abstractmethod
    def is_snv(self):
        pass

    @abstractmethod
    def normalized(self):
        pass


class BiallelicVcfVariant(VcfVariant):

    __slots__ = ("position", "reference_allele", "alternative_allele")

    def __init__(self, position: int, reference_allele: str, alternative_allele: str):
        """
        Multi-ALT sites are not modelled.
        """
        self.position = position
        self.reference_allele = reference_allele
        self.alternative_allele = alternative_allele

    def __repr__(self):
        return "BiallelicVcfVariant({}, {!r}, {!r})".format(
            self.position, self.reference_allele, self.alternative_allele
        )

    def __hash__(self):
        return hash((self.position, self.reference_allele, self.alternative_allele))

    def __eq__(self, other):
        return (
            (self.position == other.position)
            and (self.reference_allele == other.reference_allele)
            and (self.alternative_allele == other.alternative_allele)
        )

    def __lt__(self, other):
        return (self.position, self.reference_allele, self.alternative_allele) < (
            other.position,
            other.reference_allele,
            other.alternative_allele,
        )

    def get_ref_allele(self):
        return self.reference_allele

    def get_alt_allele(self):
        return self.alternative_allele

    def get_alt_allele_list(self):
        return [self.alternative_allele]

    def get_allele(self, a):
        if a == 0:
            return self.reference_allele
        elif a == 1:
            return self.alternative_allele
        else:
            raise VcfInvalidAllele(f"Querying invalid allele {a} (highest id was 1")

    def is_snv(self) -> bool:
        return (self.reference_allele != self.alternative_allele) and (
            len(self.reference_allele) == len(self.alternative_allele) == 1
        )

    def normalized(self) -> "BiallelicVcfVariant":
        """
        Return a normalized version of this variant.

        Common prefixes and/or suffixes between the reference and alternative allele are removed,
        and the position is adjusted as necessary.

        >>> BiallelicVcfVariant(100, 'GCTGTT', 'GCTAAATT').normalized()
        BiallelicVcfVariant(103, 'G', 'AAA')
        """
        pos, ref, alt = self.position, self.reference_allele, self.alternative_allele
        while len(ref) >= 1 and len(alt) >= 1 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]

        while len(ref) >= 1 and len(alt) >= 1 and ref[0] == alt[0]:
            ref, alt = ref[1:], alt[1:]
            pos += 1

        return BiallelicVcfVariant(pos, ref, alt)


class MultiallelicVcfVariant(VcfVariant):
    __slots__ = ("position", "reference_allele", "alternative_alleles")

    def __init__(self, position: int, reference_allele: str, alternative_alleles: Sequence[str]):
        self.position = position
        self.reference_allele = reference_allele
        self.alternative_alleles = tuple(alternative_alleles)

    def __repr__(self):
        return "MultiallelicVcfVariant({}, {!r}, {!r})".format(
            self.position, self.reference_allele, self.alternative_alleles
        )

    def __hash__(self):
        return hash((self.position, self.reference_allele, self.alternative_alleles))

    def __eq__(self, other):
        return (
            (self.position == other.position)
            and (self.reference_allele == other.reference_allele)
            and (len(self.alternative_alleles) == len(other.alternative_alleles))
            and all(a == b for a, b in zip(self.alternative_alleles, other.alternative_alleles))
        )

    def __lt__(self, other):
        if (self.position, self.reference_allele) != (other.position, other.reference_allele):
            return (self.position, self.reference_allele) < (other.position, other.reference_allele)
        if len(self.alternative_alleles) != len(other.alternative_alleles):
            return len(self.alternative_alleles) < len(other.alternative_alleles)
        for alt_self, alt_other in zip(
            sorted(self.alternative_alleles), sorted(other.alternative_alleles)
        ):
            if alt_self != alt_other:
                return alt_self < alt_other

        return False

    def get_ref_allele(self):
        return self.reference_allele

    def get_alt_allele(self):
        return self.alternative_alleles[0]

    def get_alt_allele_list(self):
        return self.alternative_alleles

    def get_allele(self, a):
        if a == 0:
            return self.reference_allele
        else:
            return self.alternative_alleles[a - 1]

    def is_snv(self) -> bool:
        return any(self.reference_allele != alt for alt in self.alternative_alleles) and (
            len(self.reference_allele) == 1
            and all(len(alt) == 1 for alt in self.alternative_alleles)
        )

    def normalized(self) -> "MultiallelicVcfVariant":
        """
        Return a normalized version of this variant.

        Common prefixes and/or suffixes between the reference and alternative allele are removed,
        and the position is adjusted as necessary.

        """
        pos, ref, alts = self.position, self.reference_allele, self.alternative_alleles
        while ref and all(alts) and all(ref[-1] == alt[-1] for alt in alts):
            ref, alts = ref[:-1], tuple(alt[:-1] for alt in alts)

        while ref and all(alts) and all(ref[0] == alt[0] for alt in alts):
            ref, alts = ref[1:], tuple(alt[1:] for alt in alts)
            pos += 1

        return MultiallelicVcfVariant(pos, ref, alts)


class GenotypeLikelihoods:
    __slots__ = "log_prob_genotypes"

    def __init__(self, log_prob_genotypes: List[float]):
        """Likelihoods of all genotypes to be given as log10 of
        the original probability."""
        self.log_prob_genotypes = log_prob_genotypes

    def __repr__(self):
        return f"GenotypeLikelihoods({self.log_prob_genotypes})"

    def __eq__(self, other):
        if other is None:
            return False
        if self.log_prob_genotypes is None and other.log_prob_genotypes is None:
            return True
        return self.log_prob_genotypes == other.log_prob_genotypes

    def log10_probs(self) -> List[float]:
        return self.log_prob_genotypes

    def log10_prob_of(self, genotype_index: int) -> float:
        return self.log10_probs()[genotype_index]

    def as_phred(
        self, ploidy: int = 2, regularizer: Optional[float] = None
    ) -> PhredGenotypeLikelihoods:
        if regularizer is None:
            # shift log likelihoods such that the largest one is zero
            m = max(self.log_prob_genotypes)
            return PhredGenotypeLikelihoods(
                [round((prob - m) * -10) for prob in self.log_prob_genotypes], ploidy=ploidy
            )
        else:
            p = [10**x for x in self.log_prob_genotypes]
            s = sum(p)
            p = [x / s + regularizer for x in p]
            m = max(p)
            return PhredGenotypeLikelihoods(
                [round(-10 * math.log10(x / m)) for x in p], ploidy=ploidy
            )


class VariantTable:
    """
    For a single chromosome, store variants and their genotypes.
    Each row of this table contains a variant, each column
    contains the genotypes of a single sample.

    chromosome -- chromosome name
    samples -- list of sample names
    """

    def __init__(self, chromosome: str, samples: List[str]):
        self.chromosome = chromosome
        self.samples = samples
        self.genotypes: List[List[Genotype]] = [[] for _ in samples]
        self.phases: List[List[Optional[VariantCallPhase]]] = [[] for _ in samples]
        self.allele_depths: List[List[Optional[int]]] = [[] for _ in samples]
        self.genotype_likelihoods: List[List[Optional[GenotypeLikelihoods]]] = [[] for _ in samples]
        self.variants: List[VcfVariant] = []
        self._sample_to_index = {sample: index for index, sample in enumerate(samples)}

    def __len__(self) -> int:
        return len(self.variants)

    # fmt: off
    # def add_sample(self, name, genotypes):
    # "Add a column to the table"
    # if len(genotypes) != len(self.variants):
    # raise ValueError('Expecting as many genotypes as there are variants')
    # self._name_to_index[name] = len(self.samples)
    # self.samples.append(name)
    # self.genotypes.append(genotypes)
    # fmt: on

    def add_variant(
        self,
        variant: VcfVariant,
        genotypes: Sequence[Genotype],
        phases: Sequence[Optional[VariantCallPhase]],
        genotype_likelihoods: Sequence[Optional[GenotypeLikelihoods]],
        allele_depths: Sequence[Optional[int]],
    ) -> None:
        """Add a row to the table"""
        if len(genotypes) != len(self.genotypes):
            raise ValueError("Expecting as many genotypes as there are samples")
        if len(phases) != len(self.phases):
            raise ValueError("Expecting as many phases as there are samples")
        if len(allele_depths) != len(self.allele_depths):
            raise ValueError("Expecting as many allele_depths as there are samples")
        self.variants.append(variant)
        for i, genotype in enumerate(genotypes):
            assert isinstance(genotype, Genotype)
            self.genotypes[i].append(genotype)
        for i, phase in enumerate(phases):
            self.phases[i].append(phase)
        for i, gl in enumerate(genotype_likelihoods):
            self.genotype_likelihoods[i].append(gl)
        for i, depth in enumerate(allele_depths):
            self.allele_depths[i].append(depth)

    def genotypes_of(self, sample: str) -> List[Genotype]:
        """Retrieve genotypes by sample name"""
        return self.genotypes[self._sample_to_index[sample]]

    def set_genotypes_of(self, sample: str, genotypes: List[Genotype]) -> None:
        """Set genotypes by sample name"""
        assert len(genotypes) == len(self.variants)
        self.genotypes[self._sample_to_index[sample]] = genotypes

    def genotype_likelihoods_of(self, sample: str) -> List[Optional[GenotypeLikelihoods]]:
        """Retrieve genotype likelihoods by sample name"""
        return self.genotype_likelihoods[self._sample_to_index[sample]]

    def set_genotype_likelihoods_of(
        self, sample: str, genotype_likelihoods: List[Optional[GenotypeLikelihoods]]
    ) -> None:
        """Set genotype likelihoods by sample name"""
        assert len(genotype_likelihoods) == len(self.variants)
        self.genotype_likelihoods[self._sample_to_index[sample]] = genotype_likelihoods

    def phases_of(self, sample: str) -> List[Optional[VariantCallPhase]]:
        """Retrieve phases by sample name"""
        return self.phases[self._sample_to_index[sample]]

    def num_of_blocks_of(self, sample: str) -> int:
        """Retrieve the number of blocks of the sample"""
        return len(
            {i.block_id for i in self.phases[self._sample_to_index[sample]] if i is not None}
        )

    def allele_depths_of(self, sample: str) -> List[Tuple[int, ...]]:
        """Retrieve allele depths by sample name"""
        depths: List[Tuple[int, ...]] = []
        for depth_code in self.allele_depths[self._sample_to_index[sample]]:
            assert depth_code is not None
            c = depth_code
            depth = []
            while c > 0:
                depth.append(c & 4095)
                c = c >> 12
            depths.append(tuple(depth))
        return depths

    def id_of(self, sample: str) -> int:
        """Return a unique int id of a sample given by name"""
        return self._sample_to_index[sample]

    def remove_rows_by_index(self, indices: Iterable[int]) -> None:
        """Remove variants given by their index in the variant list"""
        for i in sorted(indices, reverse=True):
            del self.variants[i]
            for gt in self.genotypes:
                del gt[i]
            for ad in self.allele_depths:
                del ad[i]
            for ph in self.phases:
                del ph[i]
            for gl in self.genotype_likelihoods:
                del gl[i]

        for gt in self.genotypes:
            assert len(self.variants) == len(gt)
        for ph in self.phases:
            assert len(self.variants) == len(ph)
        for gl in self.genotype_likelihoods:
            assert len(self.variants) == len(gl)
        assert (
            len(self.samples)
            == len(self.genotypes)
            == len(self.phases)
            == len(self.genotype_likelihoods)
        )

    def subset_rows_by_position(self, positions: Iterable[int]) -> None:
        """Keep only rows given in positions, discard the rest"""
        positions = frozenset(positions)
        to_discard = [i for i, v in enumerate(self.variants) if v.position not in positions]
        self.remove_rows_by_index(to_discard)

    def phased_blocks_as_reads(
        self,
        sample: str,
        input_variants: Iterable[VcfVariant],
        source_id: int,
        numeric_sample_id: int,
        default_quality: int = 20,
        mapq: int = 100,
        target_ploidy: int = 2,
    ):
        """
        Yields one sorted core.Read object per phased block, encoding the phase information as
        if this block was a single sequencing read. Reads are yielded in arbitrary order.

        sample -- name of sample to retrieve
        input_variants -- variants of interest, i.e. only these variants will be retrieved
        source_id -- source_id to be assigned to each read
        numeric_sample_id -- sample id to be stored in generated reads
        default_quality -- quality assigned to heterozygous with missing phasing quality
        mapq -- mapping quality for generated reads
        """
        try:
            sample_index = self._sample_to_index[sample]
        except KeyError:
            return
        input_variant_set = set(input_variants)
        read_map: Dict[int, List[Read]] = {}  # maps block_id to list of core.Read objects
        assert (
            len(self.variants)
            == len(self.genotypes[sample_index])
            == len(self.phases[sample_index])
        )
        for variant, genotype, phase in zip(
            self.variants, self.genotypes[sample_index], self.phases[sample_index]
        ):
            if len(genotype.as_vector()) != target_ploidy:
                # skip wrong ploidy
                continue
            if variant not in input_variant_set:
                continue
            if genotype.is_homozygous():
                continue
            if phase is None or phase.phase[0] is None:
                continue
            if phase.quality is None:
                quality = default_quality
            else:
                quality = phase.quality
            if phase.block_id in read_map:
                for i, allele in enumerate(phase.phase):
                    read_map[phase.block_id][i].add_variant(variant.position, allele, quality)
            else:
                read_map[phase.block_id] = []
                for i, allele in enumerate(phase.phase):
                    name = f"{sample}_phase_{i}_block_{phase.block_id}"
                    r = Read(name, mapq, source_id, numeric_sample_id)
                    r.add_variant(variant.position, allele, quality)
                    read_map[phase.block_id].append(r)
        for key, read_list in read_map.items():
            for read in read_list:
                if len(read) > 1:
                    read.sort()
                    yield read


class MixedPhasingError(Exception):
    pass


class VcfReader:
    """
    Read a VCF file chromosome by chromosome.
    """

    def __init__(
        self,
        path: Union[str, PathLike],
        only_snvs: bool = False,
        phases: bool = False,
        genotype_likelihoods: bool = False,
        ignore_genotypes: bool = False,
        ploidy: Optional[int] = None,
        mav: bool = False,
        allele_depth: bool = False,
    ):
        """
        path -- Path to VCF file
        only_snvs -- Whether to only include SNVs in the list of variants.
        ignore_genotypes -- In case of genotyping algorithm, no genotypes may be given in
                                vcf, so ignore all genotypes
        ploidy -- Ploidy of the samples
        """
        # TODO Always include deletions since they can 'overlap' other variants
        self._only_snvs = only_snvs
        self._vcf_reader = VariantFile(os.fspath(path))
        self._path = path
        self._phases = phases
        self._genotype_likelihoods = genotype_likelihoods
        self._ignore_genotypes = ignore_genotypes
        self.samples = list(self._vcf_reader.header.samples)  # intentionally public
        self.contigs = self._vcf_reader.header.contigs
        self.ploidy = ploidy
        self.mav = mav
        self.allele_depth = allele_depth
        logger.debug("Found %d sample(s) in the VCF file.", len(self.samples))

    def __enter__(self):
        return self

    def __exit__(self, *args):
        # follows same structure as for ReadSetReader
        self.close()

    def close(self):
        self._vcf_reader.close()

    @property
    def path(self) -> str:
        return self._vcf_reader.filename.decode()

    def index_exists(self) -> bool:
        """ "Check if VCF is indexed (.tbi or .csi)"""
        return self._vcf_reader.index is not None

    def _fetch(self, chromosome: str, start: int = 0, end: Optional[int] = None):
        try:
            records = self._vcf_reader.fetch(chromosome, start=start, stop=end)
        except ValueError as e:
            if "invalid contig" in e.args[0]:
                raise VcfInvalidChromosome(e.args[0]) from None
            elif "fetch requires an index" in e.args[0]:
                raise VcfIndexMissing(f"{self._path} is missing an index (.tbi or .csi)") from None
            else:
                raise
        return records

    def fetch(self, chromosome: str, start: int = 0, end: Optional[int] = None) -> VariantTable:
        """
        Fetch records from a single chromosome, optionally restricted to a single region.

        Return a VariantTable object.
        """
        records = list(self._fetch(chromosome, start=start, end=end))
        return self._process_single_chromosome(chromosome, records)

    def fetch_regions(
        self, chromosome: str, regions: Iterable[Tuple[int, Optional[int]]]
    ) -> VariantTable:
        """
        Fetch records from a single chromosome that overlap the given regions.

        :param regions: a list of start, end tuples (end can be None)
        """
        records = []
        for start, end in regions:
            records.extend(list(self._fetch(chromosome, start=start, end=end)))
        return self._process_single_chromosome(chromosome, records)

    def __iter__(self) -> Iterator[VariantTable]:
        """
        Yield VariantTable objects for each chromosome.

        Multi-ALT sites are skipped.
        """
        for chromosome, records in itertools.groupby(self._vcf_reader, lambda record: record.chrom):
            yield self._process_single_chromosome(chromosome, records)

    @staticmethod
    def _extract_HP_phase(call) -> Optional[VariantCallPhase]:
        hp = call.get("HP")
        if hp is None or hp == (".",):
            return None
        fields = [[int(x) for x in s.split("-")] for s in hp]
        for i in range(len(fields)):
            assert fields[0][0] == fields[i][0]
        block_id = fields[0][0]
        order = [field[1] - 1 for field in fields]
        phase = call["GT"]
        phase = tuple(phase[order.index(i)] for i in range(len(order)))
        return VariantCallPhase(block_id=block_id, phase=phase, quality=call.get("PQ", None))

    @staticmethod
    def _extract_GT_PS_phase(call) -> Optional[VariantCallPhase]:
        if not call.phased:
            return None
        is_het = not all(x == call["GT"][0] for x in call["GT"])
        if not is_het:
            return None
        block_id = call.get("PS", 0)
        phase = call["GT"]
        return VariantCallPhase(block_id=block_id, phase=phase, quality=call.get("PQ", None))

    @staticmethod
    def _extract_AD_depth(call) -> int:
        """
        Allele depth of sample are coded into a single integer to save space in memory.
        Encoding: Lowest 12 bits = depth of allele 0, next 12 bits = depth of allele 1, etc.
        Alleles with depth 0 are included. Maximum allele depths is 4095 (2^12 - 1).
        """
        depths = call["AD"]
        depth_code = 0
        if depths and None not in depths:
            for depth in reversed(depths):
                if depth > 4095:
                    warn_once(
                        logger, "Allele depths of 4096 or higher detected. Cutting them off to 4095"
                    )
                cnt = min(4095, depth)
                depth_code = (depth_code << 12) + cnt

        return depth_code

    def _process_single_chromosome(self, chromosome: str, records) -> VariantTable:
        phase_detected = None
        n_snvs = 0
        n_other = 0
        n_multi = 0
        table = VariantTable(chromosome, self.samples)
        prev_position = None
        for record in records:
            if not record.alts:
                continue
            if len(record.alts) > 1:
                # Multi-ALT sites are not supported, yet
                n_multi += 1
                if not self.mav or len(record.alts) >= get_max_genotype_alleles():
                    continue

            pos, ref, alts = record.start, str(record.ref), [str(alt) for alt in record.alts]
            if len(ref) == 1 and all(len(alt) == 1 for alt in alts):
                n_snvs += 1
            else:
                n_other += 1
                if self._only_snvs:
                    continue

            if (prev_position is not None) and (prev_position > pos):
                raise VcfNotSortedError(
                    "VCF not ordered: {}:{} appears before {}:{}".format(
                        chromosome, prev_position + 1, chromosome, pos + 1
                    )
                )

            if prev_position == pos:
                warn_once(
                    logger, "Skipping duplicated position %s on chromosome %r", pos + 1, chromosome
                )
                continue
            prev_position = pos

            # Read phasing information (allow GT/PS or HP phase information, but not both),
            # if requested
            if self._phases:
                phases = []
                for call in record.samples.values():
                    phase = None
                    for extract_phase, phase_name in [
                        (self._extract_HP_phase, "HP"),
                        (self._extract_GT_PS_phase, "GT_PS"),
                    ]:
                        p = extract_phase(call)
                        if p is not None:
                            if phase_detected is None:
                                phase_detected = phase_name
                            elif phase_detected != phase_name:
                                raise MixedPhasingError(
                                    "Mixed phasing information in input VCF (e.g. mixing PS "
                                    "and HP fields)"
                                )
                            phase = p
                            # check for ploidy consistency and limits
                            phase_ploidy = len(p.phase)
                            if phase_ploidy > get_max_genotype_ploidy():
                                raise PloidyError(
                                    "Ploidies higher than {} are not supported."
                                    "".format(get_max_genotype_ploidy())
                                )
                            elif p is None or p.block_id is None or p.phase is None:
                                pass
                            elif self.ploidy is None:
                                self.ploidy = phase_ploidy
                            elif phase_ploidy != self.ploidy:
                                print(f"phase= {phase}")
                                raise PloidyError(
                                    "Phasing information contains inconsistent ploidy ({} and "
                                    "{})".format(self.ploidy, phase_ploidy)
                                )
                    phases.append(phase)
            else:
                phases = [None] * len(record.samples)

            # Read genotype likelihoods, if requested
            if self._genotype_likelihoods:
                genotype_likelihoods: List[Optional[GenotypeLikelihoods]] = []
                for call in record.samples.values():
                    GL = call.get("GL", None)
                    PL = call.get("PL", None)
                    # Prefer GLs (floats) over PLs (ints) if both should be present
                    if GL is not None:
                        genotype_likelihoods.append(GenotypeLikelihoods(GL))
                    elif PL is not None:
                        likelihoods = [(pl / -10) if pl is not None else None for pl in PL]
                        genotype_likelihoods.append(GenotypeLikelihoods(likelihoods))
                    else:
                        genotype_likelihoods.append(None)
            else:
                genotype_likelihoods = [None] * len(record.samples)

            if not self._ignore_genotypes:
                # check for ploidy consistency and limits
                genotype_lists = [call.get("GT", None) for call in record.samples.values()]
                for geno in genotype_lists:
                    if geno is None or None in geno:
                        continue
                    geno_ploidy = len(geno)
                    if geno_ploidy > get_max_genotype_ploidy():
                        raise PloidyError(
                            "Ploidies higher than {} are not supported."
                            "".format(get_max_genotype_ploidy())
                        )
                    elif self.ploidy is None:
                        self.ploidy = geno_ploidy
                    elif geno_ploidy != self.ploidy:
                        raise PloidyError(
                            "Inconsistent ploidy ({} and " "{})".format(self.ploidy, geno_ploidy)
                        )

                genotypes = [genotype_code(geno_list) for geno_list in genotype_lists]
            else:
                genotypes = [Genotype([]) for _ in self.samples]
                phases = [None] * len(self.samples)

            if self.allele_depth:
                depths: List[Optional[int]] = [
                    self._extract_AD_depth(call) for call in record.samples.values()
                ]
            else:
                depths = [None] * len(record.samples)

            if len(alts) == 1:
                variant: VcfVariant = BiallelicVcfVariant(
                    position=pos, reference_allele=ref, alternative_allele=alts[0]
                )
            else:
                variant = MultiallelicVcfVariant(
                    position=pos, reference_allele=ref, alternative_alleles=alts
                )
            table.add_variant(variant, genotypes, phases, genotype_likelihoods, depths)

        logger.debug(
            "Parsed %s SNVs and %s non-SNVs. Also found %s multi-ALTs.", n_snvs, n_other, n_multi
        )

        # TODO remove overlapping variants
        return table


def remove_overlapping_calls(calls):
    """
    Filter a list of variants such that no variants overlap each other.
    This applies mainly to deletions: If they occur too close to another
    variant, the deletion and the other variant are removed.

    This function also guarantees that the positions of the returned variants
    are unique. For that, it may also remove other variants (not necessarily
    involved in a deletion).

    calls -- a list of VcfVariant objects

    Return a list of VcfVariant objects.
    """
    # TODO obviously, this is not implemented ...
    return calls


@dataclass
class VcfHeader:
    format_or_info: str
    id: str
    number: Union[str, int]
    typ: str
    description: str

    def line(self):
        return (
            "##{format_or_info}=<ID={id},Number={number},Type={typ},"
            'Description="{description}">'.format(
                format_or_info=self.format_or_info,
                id=self.id,
                number=self.number,
                typ=self.typ,
                description=self.description,
            )
        )


PREDEFINED_FORMATS = {
    "GL": VcfHeader(
        "FORMAT",
        "GL",
        "G",
        "Float",
        "Genotype Likelihood, log10-scaled likelihoods of the data given the"
        " called genotype for each possible genotype generated from the"
        " reference and alternate alleles given the sample ploidy",
    ),
    "GQ": VcfHeader("FORMAT", "GQ", 1, "Integer", "Phred-scaled genotype quality"),
    "GT": VcfHeader("FORMAT", "GT", 1, "String", "Genotype"),
    "HP": VcfHeader("FORMAT", "HP", ".", "String", "Phasing haplotype identifier"),
    "PQ": VcfHeader("FORMAT", "PQ", 1, "Float", "Phasing quality"),
    "PS": VcfHeader("FORMAT", "PS", 1, "Integer", "Phase set identifier"),
    "HS": VcfHeader("FORMAT", "HS", ".", "Integer", "Haploid phase set identifier"),
    "AD": VcfHeader("FORMAT", "AD", ".", "Integer", "Observed allele depths"),
}

PREDEFINED_INFOS = {
    "AC": VcfHeader(
        "INFO",
        "AC",
        "A",
        "Integer",
        "Allele count in genotypes, for each ALT allele, in the same order as listed",
    ),
    "AN": VcfHeader("INFO", "AN", "A", "Integer", "Total number of alleles in called genotypes"),
    "END": VcfHeader("INFO", "END", 1, "Integer", "Stop position of the interval"),
    "SVLEN": VcfHeader(
        "INFO", "SVLEN", ".", "Integer", "Difference in length between REF and ALT alleles"
    ),
    "SVTYPE": VcfHeader("INFO", "SVTYPE", 1, "String", "Type of structural variant"),
}


def augment_header(header: VariantHeader, contigs: List[str], formats: List[str], infos: List[str]):
    """
    Add contigs, formats and infos to a VariantHeader.

    formats and infos are given as a list of strings, where each item is the ID of the header
    line to add. The full header info (Number, Type, Description) is taken from the PREDEFINED_*
    constants above. Any other FORMATs or INFOs that are not predefined will raise a VcfError.

    The header is modified in place.
    """
    for contig in contigs:
        header.contigs.add(contig)

    for fmt in formats:
        if fmt in header.formats:
            header.formats[fmt].remove_header()
        try:
            h = PREDEFINED_FORMATS[fmt]
        except KeyError:
            raise VcfError(f"FORMAT {fmt!r} not defined in VCF header") from None
        header.add_line(h.line())

    for info in infos:
        try:
            h = PREDEFINED_INFOS[info]
        except KeyError:
            raise VcfError(f"INFO {info!r} not defined in VCF header") from None
        header.add_line(h.line())


def missing_headers(path: str) -> Tuple[List[str], List[str], List[str]]:
    """
    Find contigs, FORMATs and INFOs that are used within the body of a VCF file, but are
    not listed in the header or that have an incorrect type.

    Return a tuple (contigs, formats, infos) where each of the items are lists of
    strings.

    The reason this function exists is that pysam.VariantFile crashes when we
    try to write a VCF record to it that uses contigs, INFOs or FORMATs that
    are missing from the header. See also
    <https://github.com/pysam-developers/pysam/issues/771>
    """
    with VariantFile(path) as variant_file:
        header = variant_file.header.copy()
        # Check for FORMATs that do not have the expected type
        incorrect_formats = []
        for fmt, v in variant_file.header.formats.items():
            if fmt not in PREDEFINED_FORMATS:
                continue
            h = PREDEFINED_FORMATS[fmt]
            if v.number != h.number or (
                # "Float" instead of "Integer" is ok
                v.type != h.typ
                and not (v.type == "Float" and h.typ == "Integer")
            ):
                if fmt == "PS" and v.type != h.typ:
                    raise VcfError(
                        "The input VCF/BCF contains phase set ('PS') tags that are of the"
                        " non-standard type '{}' instead of 'Integer'. WhatsHap cannot"
                        " overwrite these as it could produce inconsistent files."
                        " To proceed, you can use 'whatshap unphase' to remove phasing"
                        " information from the input file".format(v.type)
                    )
                incorrect_formats.append(fmt)

        # Iterate through entire file and check which contigs, formats and
        # info fields are used
        contigs = []  # contigs encountered, in the proper order
        seen_contigs = set()
        formats = []  # FORMATs encountered, in the proper order
        seen_formats = set()
        seen_infos: Set[str] = set()  # INFOs encountered

        for record in variant_file:
            seen_infos.update(record.info)
            if record.alts is not None:
                for alt in record.alts:
                    # If there are "vague" ALT alleles such as <INS>, <DEL> etc, then
                    # the header needs to contain a LEN info entry even if LEN
                    # is never used
                    if alt.startswith("<"):
                        seen_infos.add("END")

            # For the contigs, we maintain a set *and* a list because we want to
            # keep track of the order of the contigs.
            if record.contig not in seen_contigs:
                contigs.append(record.contig)
            seen_contigs.add(record.contig)

            for fmt in record.format:
                if fmt not in seen_formats:
                    formats.append(fmt)
                seen_formats.add(fmt)

    # Determine which contigs are missing from the header
    header_contigs = set(header.contigs)
    missing_contigs = []
    for contig in contigs:
        if contig not in header_contigs:
            missing_contigs.append(contig)

    # Determine which FORMATs are missing from the header
    header_formats = set(header.formats)
    missing_formats = []
    for fmt in formats:
        if fmt in header_formats:
            continue
        missing_formats.append(fmt)

    # Determine which INFOs are missing from the header
    missing_infos = list(set(seen_infos) - set(header.info))

    return (missing_contigs, incorrect_formats + missing_formats, missing_infos)


@dataclass
class GenotypeChange:
    sample: str
    chromosome: str
    variant: VcfVariant
    old_gt: Genotype
    new_gt: Genotype


class VcfAugmenter(ABC):
    def __init__(
        self,
        in_path: str,
        command_line: Optional[str],
        out_file: TextIO = sys.stdout,
        include_haploid_phase_sets: bool = False,
    ):
        """
        in_path -- Path to input VCF, used as template.
        command_line -- A string that will be added as a VCF header entry
            (use None to not add this to the VCF header)
        out_file -- Open file-like object to which VCF is written.
        tag -- which type of tag to write, either 'PS' or 'HP'. 'PS' is standardized;
            'HP' is compatible with GATK’s ReadBackedPhasing.
        """
        # TODO This is slow because it reads in the entire VCF one extra time
        logger.debug("Reading the input VCF to find possibly missing headers")
        contigs, formats, infos = missing_headers(in_path)
        logger.debug("Missing contigs: %s", contigs)
        logger.debug("Missing formats: %s", formats)
        logger.debug("Missing infos: %s", infos)
        # TODO It would actually look nicer if the custom HS header was directly below PS
        if include_haploid_phase_sets and "HS" not in formats:
            formats.append("HS")
        # We repair the header (adding missing contigs, formats, infos) of the *input* VCF because
        # we will modify the records that we read, and these are associated with the input file.
        self._reader = VariantFile(in_path)
        augment_header(self._reader.header, contigs, formats, infos)
        if command_line is not None:
            command_line = '"' + command_line.replace('"', "") + '"'
            self._reader.header.add_meta("commandline", command_line)
        self.setup_header(self._reader.header)
        self._writer = VariantFile(out_file, mode="w", header=self._reader.header)
        self._unprocessed_record: Optional[VariantRecord] = None
        self._reader_iter = iter(self._reader)

    @abstractmethod
    def setup_header(self, header):
        pass

    def close(self):
        self._writer.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    @property
    def samples(self) -> List[str]:
        return list(self._reader.header.samples)

    def _record_modifier(self, chromosome: str):
        for record in self._iterrecords(chromosome):
            yield record
            self._writer.write(record)

    def _iterrecords(self, chromosome: str) -> Iterable[VariantRecord]:
        """Yield all records for the target chromosome"""
        n = 0
        if self._unprocessed_record is not None:
            assert self._unprocessed_record.chrom == chromosome
            yield self._unprocessed_record
            n += 1
        for record in self._reader_iter:
            n += 1
            if record.chrom != chromosome:
                # save it for later
                self._unprocessed_record = record
                assert n != 1
                return
            yield record

    def write_unchanged(self, chromosome: str) -> None:
        """
        Write all variants on one chromosome unchanged to output VCF
        chromosome -- name of chromosome
        """
        for record in self._iterrecords(chromosome):
            self._writer.write(record)


class PhasedVcfWriter(VcfAugmenter):
    """
    Read in a VCF file and write it back out with added phasing information.

    Avoid reading in full chromosomes as that uses too much memory for
    multi-sample VCFs.
    """

    def __init__(
        self,
        in_path: str,
        command_line: Optional[str],
        out_file: TextIO = sys.stdout,
        tag: str = "PS",
        ploidy: int = 2,
        include_haploid_sets: bool = False,
        only_snvs: bool = False,
        mav: bool = False,
    ):
        """
        in_path -- Path to input VCF, used as template.
        command_line -- A string that will be added as a VCF header entry
            (use None to not add this to the VCF header)
        out_file -- Open file-like object to which VCF is written.
        tag -- which type of tag to write, either 'PS' or 'HP'. 'PS' is standardized;
            'HP' is compatible with GATK’s ReadBackedPhasing.
        """
        if tag not in ("HP", "PS"):
            raise ValueError('Tag must be either "HP" or "PS"')
        self.tag = tag
        self.ploidy = ploidy
        super().__init__(in_path, command_line, out_file, include_haploid_sets)
        self._phase_tag_found_warned = False
        self._set_phasing_tags = self._set_HP if tag == "HP" else self._set_PS
        self._only_snvs = only_snvs
        self._mav = mav

    def setup_header(self, header: VariantHeader):
        """Called by baseclass constructor"""

        # FreeBayes adds phasing=none to its VCF output - remove that.
        for hr in header.records:
            if hr.key == "phasing":
                hr.remove()
                break

        header.add_line(PREDEFINED_FORMATS[self.tag].line())

    def _set_HP(
        self,
        call: VariantRecordSample,
        component: int,
        phase: Tuple[int, ...],
        haploid_component: Optional[Iterable[int]] = None,
    ):
        """
        values -- tag dict to update
        component -- name of the component
        phase -- tuple of alleles
        """
        assert all(allele in [0, 1] or self._mav for allele in phase)
        call["HP"] = ",".join(f"{component + 1}-{allele + 1}" for allele in phase)
        if haploid_component:
            call["HS"] = [comp + 1 for comp in haploid_component]

    def _set_PS(
        self,
        call: VariantRecordSample,
        component: int,
        phase: Tuple[int, ...],
        haploid_component: Optional[Iterable[int]] = None,
    ):
        """
        values -- tag dict to update
        component -- name of the component
        phase -- tuple of alleles
        """
        assert all(allele in [0, 1] or self._mav for allele in phase)
        call["PS"] = component + 1
        call["GT"] = phase
        if haploid_component:
            call["HS"] = [comp + 1 for comp in haploid_component]
        call.phased = True

    def write(
        self,
        chromosome: str,
        sample_superreads: Dict[str, ReadSet],
        sample_components: Dict,
        sample_haploid_components=None,
    ):
        """
        Add phasing information to all variants on a single chromosome.

        chromosome -- name of chromosome
        sample_superreads -- dictionary that maps a sample name to a superread
        sample_components -- a dictionary that maps a sample to its connected components

            Each component in turn is a dict that maps each variant position to a
            component, where a component is identified by the position of its
            left-most variant

        Since coordinates within the superreads are used to identify variants,
        variants at duplicate positions (allowed by the VCF spec) are currently
        not supported.

        Returns a list of changed genotyes (i.e. a list of GenotypeChange objects)
        """
        genotype_changes = []
        # TODO
        sample_phases: Dict[str, Dict] = dict()
        sample_genotypes: Dict[str, Dict] = dict()
        for sample, superreads in sample_superreads.items():
            sample_phases[sample] = {}
            sample_genotypes[sample] = {}
            for variants in zip(*superreads):
                phasing = tuple(v.allele for v in variants)
                allowed_alleles = True
                for allele in phasing:
                    if allele not in [0, 1] and not self._mav:
                        allowed_alleles = False
                        break
                if allowed_alleles:
                    sample_phases[sample][variants[0].position] = phasing
                    sample_genotypes[sample][variants[0].position] = Genotype(list(phasing))

        prev_pos = None
        for record in self._record_modifier(chromosome):
            self._remove_existing_phasing(record, list(sample_superreads))
            pos = record.start
            if not record.alts:
                continue
            if len(record.alts) > 1 and not self._mav:
                # we do not phase multiallelic sites unless requested
                continue
            if pos == prev_pos:
                # duplicate position, skip it
                continue
            is_snv = len(str(record.ref)) == 1 and len(str(record.alts[0])) == 1
            if self._only_snvs and not is_snv:
                continue

            # Determine whether the variant is phased in any sample
            for sample in self.samples:
                if sample in sample_superreads:
                    components = sample_components[sample]
                    phases = sample_phases[sample]
                    if pos in components and pos in phases:
                        break
            else:
                continue

            # Set phase tag for all target samples
            for sample in sample_superreads:
                call: VariantRecordSample = record.samples[sample]
                components = sample_components[sample]
                haploid_components = (
                    sample_haploid_components[sample] if sample_haploid_components else None
                )
                phases = sample_phases[sample]
                genotypes = sample_genotypes[sample]

                if (
                    self.tag in call
                    and call[self.tag] is not None
                    and not self._phase_tag_found_warned
                ):
                    logger.warning(
                        "Ignoring existing phasing information "
                        "found in input VCF ({} tag exists).".format(self.tag)
                    )
                    self._phase_tag_found_warned = True

                gt_type = genotype_code(call["GT"])
                is_het = not gt_type.is_homozygous()

                # is genotype to be changed?
                if pos in genotypes and genotypes[pos] != gt_type:
                    # call['GT'] = INT_TO_UNPHASED_GT[genotypes[pos]]
                    call["GT"] = tuple(genotypes[pos].as_vector())
                    variant: Union[BiallelicVcfVariant, MultiallelicVcfVariant]
                    if len(record.alts) > 1:
                        variant = MultiallelicVcfVariant(record.start, record.ref, record.alts)
                    else:
                        variant = BiallelicVcfVariant(record.start, record.ref, record.alts[0])
                    genotype_changes.append(
                        GenotypeChange(sample, chromosome, variant, gt_type, genotypes[pos])
                    )
                    is_het = not genotypes[pos].is_homozygous()

                if pos in components and pos in phases and is_het:
                    haploid_component = (
                        haploid_components[pos]
                        if (
                            haploid_components
                            and pos in haploid_components
                            and len(haploid_components[pos]) == self.ploidy
                        )
                        else None
                    )
                    self._set_phasing_tags(call, components[pos], phases[pos], haploid_component)
                else:
                    # Unphased
                    call[self.tag] = None
            prev_pos = pos
        return genotype_changes

    def _remove_existing_phasing(self, record: VariantRecord, samples: Iterable[str]):
        if self.tag == "PS":
            for sample in samples:
                call = record.samples[sample]
                if "GT" not in call:
                    continue
                call.phased = False
                if call["GT"] is not None and all(allele is not None for allele in call["GT"]):
                    call["GT"] = sorted(call["GT"])


def genotype_code(gt: Optional[Tuple[Optional[int], ...]]) -> Genotype:
    """Return genotype encoded as PyVCF-compatible number"""
    if gt is None:
        result = Genotype([])
    elif any(allele is None for allele in gt):
        result = Genotype([])
    else:
        result = Genotype([allele for allele in gt])  # type: ignore
    return result


# class to print computed genotypes,likelihoods (still needs to be improved...)
# in input vcf, currently GT is still required..


class GenotypeVcfWriter(VcfAugmenter):
    """
    Read in a VCF file and write it back out with added genotyping information.

    Avoid reading in full chromosomes as that uses too much memory for
    multi-sample VCFs.
    """

    def __init__(self, in_path: str, command_line: Optional[str], out_file: TextIO = sys.stdout):
        """
        in_path -- Path to input VCF, used as template.
        command_line -- A string that will be added as a VCF header entry.
        out_file -- Open file-like object to which VCF is written.
        """
        super().__init__(in_path, command_line, out_file)

    def setup_header(self, header: VariantHeader):
        """Called by baseclass constructor"""
        header.add_line(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype computed by WhatsHap genotyping algorithm">'
        )
        header.add_line(
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred-scaled genotype quality computed by WhatsHap genotyping algorithm">'
        )
        header.add_line(
            '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10-scaled likelihoods for genotypes: 0/0, 0/1, 1/1, computed by WhatsHap genotyping algorithm">'
        )

    def write_genotypes(
        self, chromosome: str, variant_table: VariantTable, only_snvs, ploidy: int = 2
    ) -> None:
        """
        Add genotyping information to all variants on a single chromosome.

        chromosome -- name of chromosome
        variant_table -- contains genotyping information for all accessible variant positions
        leave_unchanged -- if True, leaves records of current chromosome unchanged
        """

        # map positions to index
        genotyped_variants = dict()
        for i in range(len(variant_table)):
            genotyped_variants[variant_table.variants[i].position] = i

        # INT_TO_UNPHASED_GT = {0: (0, 0), 1: (0, 1), 2: (1, 1), -1: None}
        GT_GL_GQ = frozenset(["GT", "GL", "GQ"])
        for record in self._record_modifier(chromosome):
            pos = record.start
            if not record.alts:
                continue

            for sample, call in record.samples.items():
                geno = Genotype([])
                n_alleles = 1 + len(record.alts)
                n_genotypes = binomial_coefficient(ploidy + n_alleles - 1, n_alleles - 1)
                geno_l = [1 / n_genotypes] * int(n_genotypes)
                geno_q = None

                # for genotyped variants, get computed likelihoods/genotypes (for all others, give uniform likelihoods)
                if pos in genotyped_variants:
                    likelihoods = variant_table.genotype_likelihoods_of(sample)[
                        genotyped_variants[pos]
                    ]
                    # likelihoods can be 'None' if position was not accessible
                    if likelihoods is not None:
                        geno_l = [l for l in likelihoods]  # type: ignore
                        geno = variant_table.genotypes_of(sample)[genotyped_variants[pos]]

                # Compute GQ
                geno_index = geno.get_index()
                geno_q = sum(geno_l[i] for i in range(n_genotypes) if i != geno_index)

                # TODO default value ok?
                # store likelihoods log10-scaled

                # Temporarily overwrite the GT field with a (fake) genotype that indicates a
                # diploid sample. Otherwise, if the GT field happens to be empty, pysam
                # complains that we are setting an incorrect number of GL values.
                call["GT"] = tuple([0] * ploidy)

                call["GL"] = [max(math.log10(j), -1000) if j > 0 else -1000 for j in geno_l]
                call["GT"] = tuple(geno.as_vector())

                # store quality as phred score
                if not geno.is_none():
                    # TODO default value ok?
                    assert geno_q is not None
                    if geno_q > 0:
                        call["GQ"] = min(round(-10.0 * math.log10(geno_q)), 10000)
                    else:
                        call["GQ"] = 10000
                else:
                    call["GQ"] = None

                record.qual = None

                # delete all other genotype information that might have been present before
                for tag in set(call.keys()) - GT_GL_GQ:
                    del call[tag]
