# kate: syntax Python;
"""
Edit distance computation.

Copied from https://bitbucket.org/marcelm/sqt/src/af255d54a21815cb9a3e0b279b431a320d4626bd/sqt/_helpers.pyx
"""
from cython.view cimport array as cvarray
import cython
from libc cimport limits

@cython.boundscheck(False)
def edit_distance(s, t, int maxdiff=-1):
	"""
	Return the edit distance between the strings s and t.
	The edit distance is the sum of the numbers of insertions, deletions,
	and mismatches that is minimally necessary to transform one string
	into the other.

	If maxdiff is not -1, then a banded alignment is performed. In that case,
	the true edit distance is returned if and only if it is maxdiff or less.
	Otherwise, a value is returned that is guaranteed to be greater than
	maxdiff, but which is not necessarily the true edit distance.
	"""
	cdef int m = len(s)  # index: i
	cdef int n = len(t)  # index: j
	cdef int e = maxdiff
	cdef int i, j, start, stop, c, prev, smallest
	cdef bint match
	cdef bytes s_bytes, t_bytes
	cdef char* sv
	cdef char* tv


	# Return early if string lengths are too different
	if e != -1 and abs(m - n) > e:
		return abs(m - n)

	s_bytes = s.encode() if isinstance(s, unicode) else s
	t_bytes = t.encode() if isinstance(t, unicode) else t
	sv = s_bytes
	tv = t_bytes

	# Skip identical prefixes
	while m > 0 and n > 0 and sv[0] == tv[0]:
		sv += 1
		tv += 1
		m -= 1
		n -= 1

	# Skip identical suffixes
	while m > 0 and n > 0 and sv[m-1] == tv[n-1]:
		m -= 1
		n -= 1

	cdef int[:] costs = cvarray(shape=(m+1,), itemsize=sizeof(int), format="i")
	if e == -1:
		# Regular (unbanded) global alignment
		with nogil:
			for i in range(m + 1):
				costs[i] = i

			# compute columns of the alignment matrix (using unit costs)
			prev = 0
			for j in range(1, n+1):
				prev = costs[0]
				costs[0] += 1
				for i in range(1, m+1):
					match = sv[i-1] == tv[j-1]
					c = min(
						prev + 1 - match,
						costs[i] + 1,
						costs[i-1] + 1)
					prev = costs[i]
					costs[i] = c
	else:
		# Banded alignment
		with nogil:
			for i in range(m + 1):
				costs[i] = i
			smallest = 0
			for j in range(1, n + 1):
				stop = min(j + e + 1, m + 1)
				if j <= e:
					prev = costs[0]
					costs[0] += 1
					smallest = costs[0]
					start = 1
				else:
					start = j - e
					prev = costs[start - 1]
					smallest = maxdiff + 1
				for i in range(start, stop):
					match = sv[i-1] == tv[j-1]
					c = min(
						prev + 1 - match,
						costs[i] + 1,
						costs[i-1] + 1)
					prev = costs[i]
					costs[i] = c
					smallest = min(smallest, c)
				if smallest > maxdiff:
					break
		if smallest > maxdiff:
			return smallest
	return costs[m]

def f(l, gap_start, gap_ext):
	return gap_start + (l-1) * gap_ext

def edit_distance_affine_gap(query,ref, mismatch_cost, int gap_start=1, int gap_extend=1):
	"""
	Compute edit distance between strings s and t using affine gap costs.
	(gotoh-algorithm)

	gap_start -- cost for starting a gap
	gap_extend -- cost for extending a gap
	mismatch_cost -- list with mismatch costs for each position

	"""

	assert len(query) == len(mismatch_cost)

	cdef int m = len(query)
	cdef int n = len(ref)
	cdef int match_cost = 0
	cdef int i,j
	cdef float prev_a, prev_b, prev_c, m_c, c_a, c_b, c_c
	cdef bytes s_bytes, t_bytes
	cdef char* sv
	cdef char* tv
	cdef len_p = 0

	s_bytes = query.encode() if isinstance(query, unicode) else query
	t_bytes = ref.encode() if isinstance(ref, unicode) else ref
	sv = s_bytes
	tv = t_bytes

	# Skip identical prefixes
	while m > 0 and n > 0 and sv[0] == tv[0]:
		sv += 1
		tv += 1
		m -= 1
		n -= 1
		len_p += 1

	# Skip identical suffixes
	while m > 0 and n > 0 and sv[m-1] == tv[n-1]:
		m -= 1
		n -= 1

	# three DP tables needed
	cdef float[:] a = cvarray(shape=(m+1,), itemsize=sizeof(float), format="f")
	cdef float[:] b = cvarray(shape=(m+1,), itemsize=sizeof(float), format="f")
	cdef float[:] c = cvarray(shape=(m+1,), itemsize=sizeof(float), format="f")

	# initialize
	a[0] = 0
	b[0] = 0
	c[0] = 0

	for i in range(1,m+1):
		a[i] = limits.INT_MAX
		b[i] = f(i, gap_start, gap_extend)
		c[i] = limits.INT_MAX

	# compute tables column wise
	for j in range(1,n+1):
		prev_a = a[0]
		prev_b = b[0]
		prev_c = c[0]

		a[0] = limits.INT_MAX
		b[0] = limits.INT_MAX
		c[0] = f(j, gap_start, gap_extend)

		for i in range(1,m+1):
			m_c = mismatch_cost[i-1 + len_p]
			if sv[i-1] == tv[j-1]:
				m_c = match_cost

			c_a = min(prev_a, prev_b, prev_c) + m_c
			c_b = min(a[i-1] + gap_start, b[i-1] + gap_extend, c[i-1] + gap_start)
			c_c = min(a[i] + gap_start, b[i] + gap_start, c[i] + gap_extend)

			prev_a = a[i]
			prev_b = b[i]
			prev_c = c[i]

			a[i] = c_a
			b[i] = c_b
			c[i] = c_c
	return int(min(a[m], b[m], c[m]))
