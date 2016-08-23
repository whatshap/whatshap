from whatshap.core import Read


def verify_mec_score_and_partitioning(dp_table, reads):
	"""Confirms that the results reported by dp_table are consistent: check
	whether the reported partitioning leads to the reported MEC score."""
	superreads, transmission_vector = dp_table.get_super_reads()
	assert len(superreads) == 1
	superreads = superreads[0]
	assert len(superreads) == 2
	# create new superreads that don't contain 3s (EQUAL COST)
	new_superreads = [Read('superread0',0), Read('superread1',0)]
	assert len(superreads[0]) == len(superreads[1])
	for i in range(len(superreads[0])):
		for j in range(2):
			v = superreads[j][i]
			allele = v.allele
			if allele == 3:
				allele = j
			new_superreads[j].add_variant(v.position, allele, v.quality)
	partitioning = dp_table.get_optimal_partitioning()
	position_to_index = { variant.position: index for index, variant in enumerate(new_superreads[0]) }
	swapped = False
	mec_score = 0
	n = 0
	for read_index, read in enumerate(reads):
		cost0 = 0
		cost1 = 0
		for variant in read:
			if variant.position in position_to_index:
				if new_superreads[0][position_to_index[variant.position]].allele != variant.allele:
					cost0 = cost0 + variant.quality
				if new_superreads[1][position_to_index[variant.position]].allele != variant.allele:
					cost1 = cost1 + variant.quality
		mec_score += min(cost0, cost1)
		if cost0 == cost1:
			continue
		haplotype = 0 if (cost0 < cost1) != swapped else 1
		if partitioning[read_index] != haplotype:
			if n == 0:
				swapped = True
			else:
				assert False
		n += 1
	print('Expected MEC score: {}, obtained MEC score: {}'.format(mec_score, dp_table.get_optimal_cost()))
	assert mec_score == dp_table.get_optimal_cost()
