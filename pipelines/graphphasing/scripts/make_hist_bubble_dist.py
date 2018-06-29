import sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

in_file = sys.argv[1]

nodes = list()
with open(in_file) as fp:
	for line in fp:
		var=line.rstrip()
		nodes.append(int(var))

n, bins, patches = plt.hist(nodes, 50, facecolor='blue', alpha=0.75)

#l = plt.plot(bins)

plt.savefig("bubble_dist_his.pdf")

