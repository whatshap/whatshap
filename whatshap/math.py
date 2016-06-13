# This function was copied from Python 3.5’s statistics module.
# The StatisticsError was changed to a ValueError.
def _median(data):
	"""Return the median (middle value) of numeric data.

	When the number of data points is odd, return the middle data point.
	When the number of data points is even, the median is interpolated by
	taking the average of the two middle values:

	>>> median([1, 3, 5])
	3
	>>> median([1, 3, 5, 7])
	4.0

	"""
	data = sorted(data)
	n = len(data)
	if n == 0:
		raise ValueError("no median for empty data")
	if n%2 == 1:
		return data[n//2]
	else:
		i = n//2
		return (data[i - 1] + data[i])/2


try:
	from statistics import median
except ImportError:
	median = _median
