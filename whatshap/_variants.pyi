from typing import Iterable, Tuple, Sequence

from pysam import AlignedSegment

from .vcf import VcfVariant

def _iterate_cigar(
    variants: Sequence[VcfVariant],
    j: int,
    bam_read: AlignedSegment,
    cigartuples: Iterable[Tuple[int, int]],
) -> Iterable[Tuple[int, int, int, int]]: ...
