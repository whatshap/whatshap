import time
import logging
from collections import defaultdict
from contextlib import contextmanager
from typing import TypeVar, Iterator, Iterable, DefaultDict, Dict

logger = logging.getLogger(__name__)

T = TypeVar("T")


class StageTimer:
    """Measure run times of multiple non-overlapping stages of a program"""

    def __init__(self) -> None:
        self._start: Dict[str, float] = dict()
        self._elapsed: DefaultDict[str, float] = defaultdict(float)
        self._overall_start_time = time.time()

    def start(self, stage):
        """Start measuring elapsed time for a stage"""
        self._start[stage] = time.time()

    def stop(self, stage: str) -> float:
        """Stop measuring elapsed time for a stage."""
        t = time.time() - self._start[stage]
        if t <= 0:
            logger.warning(
                "Unreliable runtime measurements: Measured a runtime that is not positive"
            )
            t = 0
        self._elapsed[stage] += t
        del self._start[stage]
        return t

    def elapsed(self, stage: str) -> float:
        """
        Return total time spent in a stage, which is the sum of the time spans
        between calls to start() and stop(). If the timer is currently running,
        its current invocation is not counted.
        """
        return self._elapsed[stage]

    def sum(self) -> float:
        """Return sum of all times"""
        return sum(self._elapsed.values())

    def total(self) -> float:
        return time.time() - self._overall_start_time

    @contextmanager
    def __call__(self, stage: str):
        self.start(stage)
        yield
        self.stop(stage)

    def iterate(self, stage: str, iterator: Iterable[T]) -> Iterator[T]:
        """Measure iterator runtime"""

        self.start(stage)
        for item in iterator:
            self.stop(stage)
            yield item
            self.start(stage)
        self.stop(stage)
