from typing import Optional

from whatshap.core import ReadSet

def readselection(
    readset: ReadSet,
    max_cov: int,
    preferred_source_ids: Optional[set[int]] = ...,
    bridging: bool = ...,
) -> set[int]: ...
