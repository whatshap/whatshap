from abc import ABC, abstractmethod
from typing import Optional

from whatshap.core import ReadSet


class PhasingAlgorithm(ABC):
    @abstractmethod
    def get_super_reads(self) -> tuple[list[ReadSet], Optional[list[int]]]: ...

    @abstractmethod
    def get_optimal_cost(self) -> int: ...

    @abstractmethod
    def get_optimal_partitioning(self) -> list[int]: ...
