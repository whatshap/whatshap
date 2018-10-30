from itertools import permutations

def vector_error(phasing, truth, flip_cost = 1, switch_cost = 1):
	if len(phasing) != len(truth):
		raise ValueError("Incompatible phasings. Number of haplotypes is not equal ("+str(len(phasing))+" != "+str(len(truth))+").")
	num_vars = len(truth[0])
	if (num_vars == 0):
		return 0
	ploidy = len(truth)
	if (ploidy == 0):
		return 0
	for i in range(0, len(truth)):
		if len(truth[i]) != num_vars:
			raise ValueError("Inconsistent truth for phasing. Haplotypes have different lengths ( len(truth[0]="+str(num_vars)+" != len(truth["+str(i)+"]="+str(len(truth[i]))+".")
		if len(phasing[i]) != num_vars:
			raise ValueError("Inconsistent input phasing. Haplotypes have different lengths ( len(truth[0]="+str(num_vars)+" != len(phasing["+str(i)+"]="+str(len(phasing[i]))+".")
	if (len(phasing) > 6):
		print("Computing vector error with more than 6 haplotypes. This may take very long ...")

	perms = list(permutations(range(0, ploidy)))
	d = [[float("inf") for i in range(len(perms))] for j in range(num_vars)] # empty dp table
	prev = [0 for i in range(len(perms))] # auxiliary list
	prev_len = len(perms)
	
	# Initialize first column
	for i, perm in enumerate(perms):
		e = 0
		for k in range(ploidy):
			# Count flips between phasing and truth for current permutation
			e += flip_cost if (truth[k][0] != phasing[perm[k]][0] and truth[k][0] * phasing[perm[k]][0] >= 0) else 0;
		d[0][i] = e
		prev[i] = i
	
	# Iterate over all variants
	for j in range(1, num_vars):
		for i, perm in enumerate(perms):
			# Count number of flip errors if perm would be applied to this column
			flips = 0
			for k in range(ploidy):
				flips += flip_cost if (truth[k][j] != phasing[perm[k]][j] and truth[k][j] * phasing[perm[k]][j] >= 0) else 0;
				
			# Find the best previous solution by checking all rows of previous row
			min_prev_err = float("inf")
			for l, prev_perm in enumerate(perms):
				# Consider the number switches between the rows
				min_prev_err = min(min_prev_err, d[j-1][l] + switch_cost * num_switches(perm, prev_perm))
			d[j][i] = flips + min_prev_err
	
	# Vector error is smallest entry in last column
	return min(d[-1])
			
def num_switches(perm1, perm2):
	return sum([1 for i in range(len(perm1)) if perm1[i] != perm2[i]])