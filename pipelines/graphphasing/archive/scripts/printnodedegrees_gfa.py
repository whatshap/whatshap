#!/usr/bin/python

from collections import defaultdict
import fileinput
import re

degrees = defaultdict(int)

for line in fileinput.input():
	if line[0] != 'L': continue
	parts = line.split('\t')
	degrees[parts[1]] += 1
	degrees[parts[3]] += 1

print(str(sorted([(d, degrees[d]) for d in degrees], key = lambda x: -x[1])).replace('), ', ')\n').translate(None, '[](),\''))
