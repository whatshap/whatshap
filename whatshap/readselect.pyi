from whatshap.core import ReadSet

def readselection(
    readset: ReadSet,
    max_cov: int,
    preferred_source_ids: set[int] | None = ...,
    bridging: bool = ...,
) -> set[int]: ...
