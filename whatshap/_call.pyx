from libcpp.string cimport string

def enumerate_reference_kmers(string reference, int k):
	cdef int A = ord('A')
	cdef int C = ord('C')
	cdef int G = ord('G')
	cdef int T = ord('T')
	cdef c = 0
	cdef int h = 0
	cdef int mask = (1 << (2*k)) - 1
	cdef int i = 0
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
		else:
			h = ((h << 2) | 0) & mask
		if i >= k-1:
			yield (h, i+1)
