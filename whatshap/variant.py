from dataclasses import dataclass


@dataclass
class Variant:
    """A single variant on a read"""

    position: int
    allele: int
    quality: int
