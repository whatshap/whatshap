


class CovMonitor:
    #TODO: This is a most simple, naive implementation. Could do this smarter.
	#TODO : Maybe change the coverage monitor from the usage of vcf_indices to the real varaint positions

	def __init__(self, length):
		self.coverage = [0] * length

	def max_coverage_in_range(self, begin, end):
		return max(self.coverage[begin:end])

	def add_read(self, begin, end):
		for i in range(begin, end):
			self.coverage[i] += 1


