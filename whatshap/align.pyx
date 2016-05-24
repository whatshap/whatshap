# kate: syntax Python;
"""
Edit distance computation.

Copied from https://bitbucket.org/marcelm/sqt/src/af255d54a21815cb9a3e0b279b431a320d4626bd/sqt/_helpers.pyx
"""
from cython.view cimport array as cvarray
import cython


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
