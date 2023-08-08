import time
import logging
from collections import defaultdict
from contextlib import contextmanager

logger = logging.getLogger(__name__)


class StageTimer:
    """Measure run times of multiple non-overlapping stages of a program"""

    def __init__(self):
        self._start = dict()
        self._elapsed = defaultdict(float)
        self._overall_start_time = time.time()

    def start(self, stage):
        """Start measuring elapsed time for a stage"""
        self._start[stage] = time.time()

    def stop(self, stage):
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

    def elapsed(self, stage):
        """
        Return total time spent in a stage, which is the sum of the time spans
        between calls to start() and stop(). If the timer is currently running,
        its current invocation is not counted.
        """
        return self._elapsed[stage]

    def sum(self):
        """Return sum of all times"""
        return sum(self._elapsed.values())

    def total(self):
        return time.time() - self._overall_start_time

    @contextmanager
    def __call__(self, stage):
        self.start(stage)
        yield
        self.stop(stage)

    def iterate(self, stage, iterator):
        """Measure iterator runtime"""

        self.start(stage)
        for item in iterator:
            self.stop(stage)
            yield item
            self.start(stage)
        self.stop(stage)
