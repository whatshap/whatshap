from scipy.stats import binom
from math import log

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
			if (overlap[i][j-i-1] >= min_overlap):
				num_pairs += 1
				num_bases += overlap[i][j-i-1]
				avg_disagr += float(diffs[i][j-i-1])
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
	sim = [[]]
	cache = {}
	for i in range(num_reads):
		sim.append([])
		for j in range(i+1, num_reads):
			(ov, di) = (overlap[i][j-i-1], diffs[i][j-i-1])
			if (ov, di) not in cache:
				cache[(ov, di)] = logratio_sim(overlap[i][j-i-1], diffs[i][j-i-1], hammingdist_same, hammingdist_diff, min_overlap)
			sim[i].append(cache[(ov, di)])
	return sim

def partial_scoring(sim, subset):
	# sim -- a two-dimensional list with the same format as the output of the score-function
	# subset -- a list of indices indicating the subset to score
	s = sorted(subset)
	if (len(sim) <= s[-1]):
		raise ValueError("Index out of bounds: Susbet contains index "+str(s[-1])+", which is higher than the size of similarity list")
	psim = []
	for index, item1 in enumerate(s):
		row = []
		for item2 in s[index+1:]:
			row.append(sim[item1][item2-item1-1])
		psim.append(row)
	return psim

def calc_overlap_and_diffs(readset):
	num_reads = len(readset)
	overlap = [[0]*(num_reads - i - 1) for i in range(num_reads)]
	diffs = [[0]*(num_reads - i - 1) for i in range(num_reads)]
	
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
			k, l = 0, 0
			while (k < len(positions[i]) and l < len(positions[j])):
				if (positions[i][k] == positions[j][l]):
					overlap[i][j-i-1] += 1
					if (alleles[i][k] != alleles[j][l]):
						diffs[i][j-i-1] += 1
					k += 1
					l += 1
				elif (positions[i][k] < positions[j][l]):
					k += 1
				else:
					l += 1
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
						total_overlap += overlap[k][l-k-1]
						total_diff += diffs[k][l-k-1]
					elif (k > l):
						total_overlap += overlap[l][k-l-1]
						total_diff += diffs[l][k-l-1]
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