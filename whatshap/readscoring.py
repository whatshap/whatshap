from scipy.stats import binom

def score(readset, ploidy, errorrate, min_overlap):
	num_reads = len(readset)
	#num_vars = readset.positions()
	num_vars = 0
	m = []

	# Map variant positions to [0,l)
	index = {}
	for position in readset.get_positions():
		index[position] = num_vars
		num_vars += 1

	begins = [num_vars+1]*num_reads
	ends = [-1]*num_reads

	# Translate to matrix and determine start and end position (both inclusive)
	for i in range(num_reads):
		line = [-1]*num_vars
		for variant in readset[i]:
			if (variant.position in index):
				begins[i] = min(begins[i], index[variant.position])
				ends[i] = max(ends[i], index[variant.position])
				line[index[variant.position]] = variant.allele
		m.append(line)

	# Calculate overlap and matches. Overlap for reads i and j is saved in overlap[min(i,j)][max(i,j)-min(i,j)-1].
	overlap = [[]]
	matches = [[]]
	for i in range(num_reads):
		overlap.append([])
		matches.append([])
		for j in range(i+1, num_reads):
			overlap[i].append(0)
			matches[i].append(0)
			start = max(begins[i], begins[j])
			end = max(ends[i], ends[j])
			overlap[i][j-i-1] = end - start + 1
			for k in range(start, end+1):
				if (m[i][k] == m[j][k]):
					matches[i][j-i-1] += 1

	del m

	# Estimate hamming distance between read pairs of same/different haplotype
	avg_disagr = 0.0
	num_pairs = 0
	num_bases = 0
	hammingdist_same = 2*(1.0-errorrate)*errorrate

	for i in range(num_reads):
		for j in range(i+1, num_reads):
			if (overlap[i][j-i-1] >= 1):
				num_pairs += 1
				num_bases += overlap[i][j-i-1]
				avg_disagr += float(overlap[i][j-i-1]-matches[i][j-i-1])

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
			sim[i].append(logratio_sim(overlap[i][j-i-1], matches[i][j-i-1], hammingdist_same, hammingdist_diff, min_overlap))

	return sim

def logratio_sim(overlap, matches, dist_same, dist_diff, min_overlap):
	if (overlap < min_overlap):
		return 0

	p_same = binom.pmf(overlap, overlap - matches, dist_same)
	p_diff = binom.pmf(overlap, overlap - matches, dist_diff)
	score = 0.0
	if (p_same == 0):
		score = -float("inf")
	elif (p_diff == 0):
		score = float("inf")
	else:
		score = 1.0
		#score = math.log(p_same / p_diff)
	
	return score

if __name__ == "__main__":
	main()
