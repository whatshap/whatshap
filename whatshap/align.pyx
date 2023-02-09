# cython: language_level=3

"""
Edit distance computation.

Copied from https://bitbucket.org/marcelm/sqt/src/af255d54a21815cb9a3e0b279b431a320d4626bd/sqt/_helpers.pyx
"""
from cython.view cimport array as cvarray
import cython
from libc cimport limits
from libc.math cimport log10
from libcpp.string cimport string
import collections


@cython.boundscheck(False)
def match_score(str alpha, str beta, float mismatch_penalty):
	cdef float matching = 0
	if alpha == beta:
		return matching
	else:
		return mismatch_penalty
def needle(seq1, seq2, dictionary, int gap_penalty, int k):
	cdef int i,j,x
	cdef int m = len(seq1)
	cdef int n = len(seq2)
	cdef int[:,:] score
	cdef float mismatching
	if seq1==seq2:

		return 0
	else:	
		x=0 
		# Skip identical prefixes
		while x < m and x < n and seq1[x] == seq2[x]:
			x+=1       #recording where the suffix match stopped

		# Skip identical suffixes
		while m > x and n > x and seq1[m-1] == seq2[n-1]:
			m -= 1
			n -= 1
		#now remove the suffix indices that we already have seen
		m-=x 
		n-=x
		score = cvarray(shape=(m+1,n+1), itemsize=sizeof(int), format="i")

		for i in range(0, m + 1):
			score[i][0] = gap_penalty * i
		for j in range(0, n + 1):
			score[0][j] = gap_penalty *j
		for i in range(1, m + 1):
			for j in range (1,n+1):
				if seq1[i-1+x]==seq2[j-1+x]:
					match = score[i - 1][j - 1] #no penalty
				else:
					if (seq1[i-1+x],seq2[j-1+x]) in dictionary:
						mismatching = float(dictionary[(seq1[i-1+x],seq2[j-1+x])])
						#print(seq1[i-1+x], '\t',seq2[j-1+x], '\t',mismatching)
					elif (seq1[i-1+x],-5) in dictionary:

						mismatching= float(dictionary[(seq1[i-1+x],-5)])
						#print(seq1[i-1+x], '\t','epsilon', '\t',mismatching)
					else:
						mismatching= float('inf')
						#print(seq1[i-1+x], '\t',seq2[j-1+x], '\t',mismatching)
					match = score[i - 1][j - 1] + mismatching

				delete = score[i - 1][j] + gap_penalty
				insert = score[i][j - 1] + gap_penalty
				score[i][j] = min(match, delete, insert)
		return score[m][n]




def needle_temp(sv, tv, dictionary, int gap_penalty, int k):

	cdef int m = len(sv)  # index: i
	cdef int n = len(tv)  # index: j
	cdef bint match
	cdef int i,j, x
	cdef int[:] costs
	if sv==tv:

		return 0
	else:	
		x=0 
		# Skip identical prefixes
		while x < m and x < n and sv[x] == tv[x]:
			x+=1       #recording where the suffix match stopped

		# Skip identical suffixes
		while m > x and n > x and sv[m-1] == tv[n-1]:
			m -= 1
			n -= 1
		#now remove the suffix indices that we already have seen
		m-=x 
		n-=x
		costs = cvarray(shape=(m+1,), itemsize=sizeof(int), format="i")
		# Regular (unbanded) global alignment
		for i in range(m + 1):
			costs[i] = i*gap_penalty

		# compute columns of the alignment matrix (using unit costs)
		prev = 0
		for j in range(1, n+1):
			prev = costs[0]
			costs[0] += gap_penalty
			for i in range(1, m+1):
				if (sv[i-1+x], tv[j-1+x]) in dictionary:
					mismatching = float(dictionary[(sv[i-1+x], tv[j-1+x])])

				elif (sv[i-1+x],-5) in dictionary:
					mismatching= float(dictionary[(sv[i-1+x],-5)])

				else:
					mismatching= float('inf')

				if sv[i-1+x] == tv[j-1+x]:
					match = 0 #no penalty
				else:
					match = prev + mismatching
					c = min(
						prev + match,
						costs[i] + gap_penalty,
						costs[i-1] + gap_penalty)
					prev = costs[i]
					costs[i] = c
		return costs[m]
	
def split(sequence, int k):
	cdef:
		kmer_list= []
		int i=0
		shortstring= ""
	if len(sequence) <= k:
		kmer_list.append(sequence)
		return kmer_list
	else:
		shortstring=sequence[i:i+k]
		kmer_list.append(shortstring)
		i+=1
		while i<=len(sequence)-k:
			shortstring=shortstring[1:]+sequence[i+k-1]
			kmer_list.append(shortstring)
			i+=1
		return kmer_list

def enumerate_all_kmers(string reference, int k):
	cdef int A = ord('A')
	cdef int C = ord('C')
	cdef int G = ord('G')
	cdef int T = ord('T')
	cdef c = 0
	cdef int h = 0
	cdef int mask = (1 << (2*k)) - 1
	cdef int i = 0
	cdef kmer_list= collections.deque([])
	for i in range(len(reference)):
		c = reference[i]
		if c == A:
			h = ((h << 2) | 0) & mask
		elif c == C:
			h = ((h << 2) | 1) & mask
		elif c == G:
			h = ((h << 2) | 2) & mask
		elif c == T:
			h = ((h << 2) | 3) & mask
		if i >= k-1:
			kmer_list.append(h)

	return kmer_list
	
def enumerate_a_kmer(str kmer, int k):
	cdef c = 0
	cdef int h = 0
	cdef int mask = (1 << (2*k)) - 1
	cdef int i = 0
	if kmer=='epsilon':
		h= -5
	else:
		for i in range(len(kmer)):
			c = kmer[i]
			if c == 'A':
				h = ((h << 2) | 0) & mask	
			elif c == 'C':
				h = ((h << 2) | 1) & mask
			elif c == 'G':
				h = ((h << 2) | 2) & mask
			elif c == 'T':
				h = ((h << 2) | 3) & mask
	return h
                
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
