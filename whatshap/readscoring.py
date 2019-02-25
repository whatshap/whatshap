from scipy.stats import binom
from math import log, floor, ceil
from .core import ReadScoring, TriangleSparseMatrix
import itertools as it
import random as rand
import logging

logger = logging.getLogger(__name__)

class SparseTriangleMatrix:
	__slots__ = ('size', 'values')
	
	def __init__(self):
		self.size = 0
		self.values = {}

	def __hash__(self):
		return hash((self.size, self.values))

	def __iter__(self):
		for i in range(self.size):
			for j in range(i):
				if self.get(i, j) != 0:
					yield (i, j, self.get(i, j))
	
	def __eq__(self, other):
		return (self.size == other.size) and \
		       (cmp(self.values, other.values) == 0)

	def __lt__(self, other):
		return (self.size, self.values) < (other.size, other.values)
	
	def entry_index(self, i, j):
		if (i < j):
			return self.entry_index(j, i)
		elif (i > j):
			return (i*(i-1)/2)+j
		else:
			return -1
		
	def get_size(self):
		return self.size

	def get(self, i, j):
		index = self.entry_index(i, j)
		assert index >= 0
		if index in self.values:
			return self.values[index]
		else:
			return 0
		
	def set(self, i, j, value):
		index = self.entry_index(i, j)
		assert index >= 0
		if index >= 0:
			self.values[index] = value
			self.size = max(i, self.size)
			self.size = max(j, self.size)

def score_global(readset, ploidy, min_overlap):
	read_scoring = ReadScoring()
	return read_scoring.scoreReadsetGlobal(readset, min_overlap, ploidy)
	
def score_local(readset, ploidy, min_overlap):
	read_scoring = ReadScoring()
	return read_scoring.scoreReadsetLocal(readset, min_overlap, ploidy)

def score_local_patternbased(readset, ploidy, errorrate, min_overlap, window_size):
	read_scoring = ReadScoring()
	return read_scoring.scoreReadsetPatterns(readset, min_overlap, ploidy, errorrate, window_size)

def partial_scoring(sim, subset):
	# sim -- a two-dimensional list with the same format as the output of the score-function
	# subset -- a list of indices indicating the subset to score
	s = sorted(subset)
	if (sim.get_size() < s[-1]):
		raise ValueError("Index out of bounds: Susbet contains index "+str(s[-1])+", which is higher than the size of similarity matrix")
		
	psim = SparseTriangleMatrix()	
	for i in range(len(s)):
		for j in range(i+1, len(s)):
			psim.set(i, j, sim.get(s[i], s[j]))
	return psim

def calc_overlap_and_diffs(readset):
	num_reads = len(readset)
	overlap = SparseTriangleMatrix()
	diffs = SparseTriangleMatrix()
	
	# Copy information from readset into lists, because direct access is very slow
	begins = [readset[i][0].position for i in range(num_reads)]
	ends = [readset[i][-1].position for i in range(num_reads)]
	positions = [[read[k].position for k in range(len(read))] for read in readset]
	alleles = [[read[k].allele for k in range(len(read))] for read in readset]
	
	# Iterate over read pairs
	for i in range(num_reads):
		for j in range(i+1, num_reads):
			# if reads do not overlap, leave overlap and differences at 0 and skip
			if (ends[i] < begins[j] or ends[j] < begins[i]):
				continue

			# perform a zigzag search over the variants of both reads
			ov = 0
			di = 0
			k, l = 0, 0
			while (k < len(positions[i]) and l < len(positions[j])):
				if (positions[i][k] == positions[j][l]):
					ov += 1
					if (alleles[i][k] != alleles[j][l]):
						di += 1
					k += 1
					l += 1
				elif (positions[i][k] < positions[j][l]):
					k += 1
				else:
					l += 1
			overlap.set(i, j, ov)
			diffs.set(i, j, di)
	return overlap, diffs

def logratio_sim(overlap, diffs, dist_same, dist_diff, min_overlap):
	if (overlap < min_overlap):
		return 0

	p_same = binom.pmf(diffs, overlap, dist_same)
	p_diff = binom.pmf(diffs, overlap, dist_diff)
	score = 0.0
	if (p_same == 0):
		score = -float("inf")
	elif (p_diff == 0):
		score = float("inf")
	else:
		score = log(p_same / p_diff)
	return score

def print_dist_between_haps(readset, overlap, diffs):
	haps = [[j for j in range(len(readset)) if parse_haplotype(readset[j].name) == i] for i in range(0, 4)]
	for i in range(4):
		for j in range(i, 4):
			total_overlap = 0
			total_diff = 0
			for k in haps[i]:
				for l in haps[j]:
					if (k < l):
						total_overlap += overlap.get(k, l)
						total_diff += diffs.get(k, l)
					elif (k > l):
						total_overlap += overlap.get(l, k)
						total_diff += diffs.get(l, k)
			print("Diff "+str(i)+" vs "+str(j)+" = "+str(total_diff/total_overlap))
			
def parse_haplotype(name):
	tokens = name.split("_")
	if (tokens[-2] == "HG00514" and tokens[-1] == "HAP1"):
		return 0
	elif (tokens[-2] == "HG00514" and tokens[-1] == "HAP2"):
		return 1
	elif (tokens[-2] == "NA19240" and tokens[-1] == "HAP1"):
		return 2
	elif (tokens[-2] == "NA19240" and tokens[-1] == "HAP2"):
		return 3
	return -1