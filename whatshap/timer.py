import time
from collections import defaultdict
from contextlib import contextmanager


class StageTimer:
	"""Measure run times of different stages of the program"""

	def __init__(self):
		self._start = dict()
		self._elapsed = defaultdict(float)

	def start(self, stage):
		"""Start measuring elapsed time for a stage"""
		self._start[stage] = time.time()

	def stop(self, stage):
		"""Stop measuring elapsed time for a stage."""
		t = time.time() - self._start[stage]
		self._elapsed[stage] += t
		return t

	def elapsed(self, stage):
		"""
		Return total time spent in a stage, which is the sum of the time spans
		between calls to start() and stop(). If the timer is currently running,
		its current invocation is not counted.
		"""
		return self._elapsed[stage]

	def total(self):
		"""Return sum of all times"""
		return sum(self._elapsed.values())

	@contextmanager
	def __call__(self, stage):
		self.start(stage)
		yield
		self.stop(stage)
