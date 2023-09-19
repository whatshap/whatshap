from typing import Iterable, Tuple, Sequence, List

from pysam import AlignedSegment

from .variants import VariantProgress
from .vcf import VcfVariant

def _iterate_cigar(
    variants: Sequence[VcfVariant],
    j: int,
    bam_read: AlignedSegment,
    cigartuples: Iterable[Tuple[int, int]],
) -> Iterable[Tuple[int, int, int, int]]: ...
def _detect_alleles(
    variants: List[VcfVariant],
    var_progress: List[VariantProgress],
    first: int,
    bam_read: AlignedSegment,
) -> Iterable[Tuple[int, str, float]]: ...
