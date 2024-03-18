from abc import ABC, abstractmethod
from typing import Tuple, List, Optional

from whatshap.core import ReadSet


class PhasingAlgorithm(ABC):
    @abstractmethod
    def get_super_reads(self) -> Tuple[List[ReadSet], Optional[List[int]]]: ...

    @abstractmethod
    def get_optimal_cost(self) -> int: ...

    @abstractmethod
    def get_optimal_partitioning(self) -> List[int]: ...
