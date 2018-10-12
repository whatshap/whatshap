from scipy.stats import binom
from math import log

def score(readset, ploidy, errorrate, min_overlap):
	num_reads = len(readset)

	# Calculate overlap and differences. Overlap for reads i and j is saved in overlap[min(i,j)][max(i,j)-min(i,j)-1].
	overlap, diffs = calc_overlap_and_diffs(readset)

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

	avg_disagr = avg_disagr / num_bases
		
	frac_same = 0.0
	if (ploidy > 1):
		frac_same = 1.0 / float(ploidy) # probability that a random read pair is from same haplotype
		
	# (frac_same*num_pairs)*errorrate + ((1-frac_same)*num_pairs)*x = num_pairs*avg_disagr
	# => ((1-frac_same)*num_pairs)*x = num_pairs*avg_disagr - (frac_same*num_pairs)*errorrate
	# => x = (num_pairs*avg_disagr - (frac_same*num_pairs)*errorrate) / ((1-frac_same)*num_pairs)
	hammingdist_diff = (num_pairs*avg_disagr - (frac_same*num_pairs)*hammingdist_same) / ((1.0-frac_same)*num_pairs)

	# Calculate the actual similarities
	sim = [[]]
	for i in range(num_reads):
		sim.append([])
		for j in range(i+1, num_reads):
			sim[i].append(logratio_sim(overlap[i][j-i-1], diffs[i][j-i-1], hammingdist_same, hammingdist_diff, min_overlap))
	return sim

def calc_overlap_and_diffs(readset):
	num_reads = len(readset)
	overlap = [[]]
	diffs = [[]]
	for i in range(num_reads):
		overlap.append([])
		diffs.append([])
		for j in range(i+1, num_reads):
			overlap[i].append(0)
			diffs[i].append(0)
			# if reads do not overlap, leave overlap and differences at 0 and skip
			if (readset[i][-1].position < readset[j][0].position or readset[j][-1].position < readset[i][0].position):
				continue

			# perform a zigzag search over the variants of both reads
			k, l = 0, 0
			while (k < len(readset[i]) and l < len(readset[j])):
				if (readset[i][k].position == readset[j][l].position):
					overlap[i][j-i-1] += 1
					if (readset[i][k].allele != readset[j][l].allele):
						diffs[i][j-i-1] += 1
					k += 1
					l += 1
				elif (readset[i][k].position < readset[j][l].position):
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
