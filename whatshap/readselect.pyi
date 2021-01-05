from typing import Optional, Set

from whatshap.core import ReadSet

def readselection(
    readset: ReadSet,
    max_cov: int,
    preferred_source_ids: Optional[Set[int]] = ...,
    bridging: bool = ...,
) -> Set[int]: ...
