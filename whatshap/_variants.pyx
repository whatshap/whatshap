def _iterate_cigar(variants, int j, bam_read, cigartuples):
	"""
	Iterate over the CIGAR of the given bam_read and variants[j:] in lockstep.

	Yield tuples (index, i, consumed, query_pos) where index is into the variants list

	i and consumed describe the split position in the cigar

	bam_read -- a pysam.AlignedSegment
	variants -- list of variants (VcfVariant objects)
	j -- index of the first variant (in the variants list) to check
	"""
	cdef:
		int ref_pos = bam_read.reference_start  # position relative to reference
		int query_pos = 0  # position relative to read
		int cigar_op
		int length
		int i
		int n = len(variants)
		int v_position

	# Skip variants that are located to the left of the read
	while j < n and variants[j].position < ref_pos:
		j += 1

	# Iterate over the CIGAR sequence (defining the alignment) and variant list in lockstep
	for i, (cigar_op, length) in enumerate(cigartuples):
		# The mapping of CIGAR operators to numbers is:
		# MIDNSHPX= => 012345678
		if j < n:
			v_position = variants[j].position
		if cigar_op in (0, 7, 8):  # M, X, = operators (match)
			# Iterate over all variants that are in this matching region
			while j < n and v_position < ref_pos + length:
				assert v_position >= ref_pos
				yield (j, i, v_position - ref_pos, query_pos + v_position - ref_pos)
				j += 1
				if j < n:
					v_position = variants[j].position
			query_pos += length
			ref_pos += length
		elif cigar_op == 1:  # I operator (insertion)
			# TODO it should work to *not* handle the variant here, but at the next M or D region
			if j < n and v_position == ref_pos:
				yield (j, i, 0, query_pos)
				j += 1
				if j < n:
					v_position = variants[j].position
			query_pos += length
		elif cigar_op == 2:  # D operator (deletion)
			# Iterate over all variants that are in this deleted region
			while j < n and v_position < ref_pos + length:
				assert v_position >= ref_pos
				yield (j, i, v_position - ref_pos, query_pos)
				j += 1
				if j < n:
					v_position = variants[j].position
			ref_pos += length
		elif cigar_op == 3:  # N operator (reference skip)
			# Iterate over all variants that are in this skipped region
			while j < n and v_position < ref_pos + length:
				assert v_position >= ref_pos
				j += 1
				if j < n:
					v_position = variants[j].position
			ref_pos += length
		elif cigar_op == 4:  # S operator (soft clipping)
			query_pos += length
		elif cigar_op == 5 or cigar_op == 6:  # H or P (hard clipping or padding)
			pass
		else:
			raise ValueError("Unsupported CIGAR operation: {}".format(cigar_op))
