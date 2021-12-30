from dataclasses import dataclass
from typing import Iterable


@dataclass
class Variant:
    """A single variant on a read"""

    position: int
    allele: int
    emission: Iterable[float]
    quality: int
