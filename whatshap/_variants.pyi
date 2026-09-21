from collections.abc import Iterable, Sequence

from pysam import AlignedSegment

from .variants import VariantProgress
from .vcf import VcfVariant

def _iterate_cigar(
    variants: Sequence[VcfVariant],
    j: int,
    bam_read: AlignedSegment,
    cigartuples: Iterable[tuple[int, int]],
) -> Iterable[tuple[int, int, int, int]]: ...
def _detect_alleles(
    variants: list[VcfVariant],
    var_progress: list[VariantProgress],
    first: int,
    bam_read: AlignedSegment,
) -> Iterable[tuple[int, str, float]]: ...
