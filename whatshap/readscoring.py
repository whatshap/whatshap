from scipy.stats import binom
from math import log, floor, ceil
from .core import ReadScoring, TriangleSparseMatrix
import itertools as it
import random as rand

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

def score(readset, ploidy, errorrate, min_overlap):
	
	read_scoring = ReadScoring()
	
	return read_scoring.scoreReadset(readset, errorrate, min_overlap, ploidy)
	
def unused(readset, ploidy, errorrate, min_overlap):
	
	num_reads = len(readset)

	# Calculate overlap and differences. Overlap for reads i and j is saved in overlap[min(i,j)][max(i,j)-min(i,j)-1].
	overlap, diffs = calc_overlap_and_diffs(readset)
	
	# Print actual hamming distance between haps
	#print_dist_between_haps(readset, overlap, diffs)

	# Estimate hamming distance between read pairs of same/different haplotype
	avg_disagr = 0.0
	num_pairs = 0
	num_bases = 0
	hammingdist_same = 2*(1.0-errorrate)*errorrate

	for i in range(num_reads):
		for j in range(i+1, num_reads):
			if (overlap.get(i, j) >= min_overlap):
				num_pairs += 1
				num_bases += overlap.get(i, j)
				avg_disagr += float(diffs.get(i, j))
	if (ploidy > 1):
		frac_same = 1.0 / float(ploidy) # probability that a random read pair is from same haplotype
	else:
		frac_same = 0.0

	if (num_bases == 0):
		avg_disagr = 0.0
		hammingdist_diff = (1.0+hammingdist_same)/2
	else:
		avg_disagr = avg_disagr / num_bases
		# (frac_same*num_pairs)*errorrate + ((1-frac_same)*num_pairs)*x = num_pairs*avg_disagr
		# => ((1-frac_same)*num_pairs)*x = num_pairs*avg_disagr - (frac_same*num_pairs)*errorrate
		# => x = (num_pairs*avg_disagr - (frac_same*num_pairs)*errorrate) / ((1-frac_same)*num_pairs)
		hammingdist_diff = (num_pairs*avg_disagr - (frac_same*num_pairs)*hammingdist_same) / ((1.0-frac_same)*num_pairs)
		hammingdist_diff = max(hammingdist_same, min((1.0+hammingdist_same)/2, hammingdist_diff))
	
	# Calculate the actual similarities
	sim = SparseTriangleMatrix()
	cache = {}
	for i in range(num_reads):
		for j in range(i+1, num_reads):
			(ov, di) = (overlap.get(i, j), diffs.get(i, j))
			if (ov < min_overlap):
				continue
			if (ov, di) not in cache:
				cache[(ov, di)] = logratio_sim(ov, di, hammingdist_same, hammingdist_diff, min_overlap)
			sim.set(i, j, cache[(ov, di)])
	return sim

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

def locality_sensitive_score(readset, ploidy, min_overlap):
	read_scoring = ReadScoring()
	return read_scoring.scoreReadsetLocally(readset, min_overlap, ploidy)
	
def also_unused(readset, ploidy, min_overlap):
	num_reads = len(readset)
	
	overlap, diffs = calc_overlap_and_diffs(readset)
	
	# Copy information from readset into lists, because direct access is very slow
	begins = [readset[i][0].position for i in range(num_reads)]
	ends = [readset[i][-1].position for i in range(num_reads)]
	positions = [[read[k].position for k in range(len(read))] for read in readset]
	alleles = [[read[k].allele for k in range(len(read))] for read in readset]
	
	# For each read r, store all reads s, which overlap with r on at least min_overlap positions
	partners = [set() for i in range(num_reads)]
	for i in range(num_reads - 1):
		for j in range(i+1, num_reads):
			if overlap.get(i, j) >= min_overlap:
				partners[i].add(j)
				partners[j].add(i)
	
	# Sparse matrix to store similarities
	sim = SparseTriangleMatrix()
	
	# Loop over all overlapping read pairs
	for i in range(num_reads - 1):
		#print("Scoring read "+str(i))
		for j in partners[i]:
			if j <= i:
				continue
				
			min_pos = max(begins[i], begins[j])
			max_pos = min(ends[i], ends[j])
			
			# Build set of all involved reads, i.e. all reads that overlap with both i and j
			co_overlappers = partners[i].intersection(partners[j])
			#print("Pair "+str(i)+"+"+str(j)+": "+str(co_overlappers))
			co_overlappers.add(i)
			co_overlappers.add(j)
			
			# Compute all relative differences between all combinations of involved reads
			relative_diffs = []
			for k, l in it.combinations(co_overlappers, 2):
				
				# if reads do not overlap each other and the target region, skip
				if overlap.get(k, l) < min_overlap or ends[k] <= min_pos or ends[l] <= min_pos or begins[k] >= max_pos or begins[l] >= max_pos:
					continue

				# perform a zigzag search over the variants of both reads
				ov = 0
				di = 0
				r, s = 0, 0
				while positions[k][r] < min_pos:
					r += 1
				while positions[l][s] < min_pos:
					s += 1
				
				while r < len(positions[k]) and s < len(positions[l]) and positions[k][r] <= max_pos and positions[l][s] <= max_pos:
					if (positions[k][r] == positions[l][s]):
						ov += 1
						if (alleles[k][r] != alleles[l][s]):
							di += 1
						r += 1
						s += 1
					elif (positions[k][r] < positions[l][s]):
						r += 1
					else:
						s += 1

				if ov >= min_overlap:
					relative_diffs.append(di / ov)
								
				#if overlap.get(k, l) >= min_overlap:
				#	relative_diffs.append(diffs.get(k, l) / overlap.get(k, l))
				
			# The bottom (1/ploidy)-fraction of differences is assumed to correspond to the (1/ploidy)-fraction of read pairs, originating from the same haplotype. We take the middle out of this.
			# The rest is divided among all other combination of haplotypes. These are (ploidy-1)*ploidy/2 many. We want to compare the relative difference of read pairs from equal haplotypes to the
			# relative difference of read pairs from the pair of haplotypes, which differ the least.
			relative_diffs.sort()
			
			#print("Relative diffs = "+str(relative_diffs))
			quantile_same = (1 / ploidy)
			quantile_diff = (1 / ploidy) + (2 / ((ploidy-1)*ploidy)) * (1 - 1 / ploidy)
			#print("Quantiles: "+str(quantile_same)+" and "+str(quantile_diff))
			quantile_same = ceil(quantile_same * (len(relative_diffs) - 1))
			quantile_diff = ceil(quantile_diff * (len(relative_diffs) - 1))
			if quantile_diff > quantile_same:
				mean_same = sum(relative_diffs[:quantile_same]) / quantile_same
				mean_diff = sum(relative_diffs[quantile_same:quantile_diff]) / (quantile_diff - quantile_same)
				#print("Quantiles: "+str(mean_same)+" and "+str(mean_diff))

				# Only compute sim, if same and different haplotype are actually distinguishable
				if mean_same < mean_diff:
					sim.set(i, j, 1 if abs(diffs.get(i, j) / overlap.get(i, j) - mean_same) < abs(diffs.get(i, j) / overlap.get(i, j) - mean_diff) else -1)
					#print("Sim: "+str(sim.get(i, j))+" ("+str(diffs.get(i, j))+" / "+str(overlap.get(i, j))+")")
				
	return sim

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