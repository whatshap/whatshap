from scipy.stats import binom
from math import log

class SparseTriangleMatrix:
	__slots__ = ('size', 'values')
	
	def __init__(self):
		self.size = 0
		self.values = {}

	def __hash__(self):
		return hash((self.size, self.values))

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