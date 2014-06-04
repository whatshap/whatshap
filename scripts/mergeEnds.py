#!/usr/bin/env python
from __future__ import print_function

# quick script to merge the (sorted) ends output to get a reads output
# -- Murray Patterson

"""
Marcel's comments:
- snps.keys() orders by string value, is that intended?

"""

import sys


def parse_line(line):
	t = line.split()
	name = t[0]
	count = 0
	mapq = 0
	is_unique = 0
	snps = {}
	for i in range(len(t)):
		if t[i] == ":":
			snps[t[i+1]] = [t[i+2], t[i+3], t[i+4]]
		if t[i] == "#":
			count = int(t[i+1])
			mapq = int(t[i+2])  # mapping quality
			is_unique = t[i+3]  # unique flag ('U' is for unique, 'R' is for repetitive, adopted from BWA XT tag)
	assert count == len(snps)
	return name, count, mapq, is_unique, snps


def main():
	if len(sys.argv) < 2 :
		print("usage : " + str(sys.argv[0]) + " endsFile")
		sys.exit(0)

	f = open(sys.argv[1],"r")

	# get first end
	e = f.readline()
	if not e :
		print("file is empty")
		sys.exit(0)

	# snps maps a position to a list [base, allele, quality]
	# parse the first line
	name, count, mapq, is_unique, snps = parse_line(e)

	# parse the remaining lines
	while True:
		ep = f.readline() # get second end
		if not ep: # no ep: end e is unpaired
			# seems we are at EOF
			if count > 1: # so simply print end e
				for p in sorted(snps.keys()) : # careful: the default is lists in no particular order, but we want snps to be ordered on their fragment
					print(p + " " + snps[p][0] + " " + snps[p][1] + " " + snps[p][2] + " : ", end='')
				print("# " + str(mapq) + " : " + is_unique)
			break

		# everything with the 'p' suffix is from the second (paired) read
		np, cp, mp, up, sp = parse_line(ep)

		if name == np: # end e pairs up with end ep
			# output merged pair (a read)
			uup = 0
			if count + cp > 1:
				for p in sorted(snps.keys()):
					print(p + " " + snps[p][0] + " " + snps[p][1] + " " + snps[p][2] + " : ", end='')
				print("-- : ", end='') # add a symbol for gap in paired-end reads
				for p in sorted(sp.keys()) :
					print(p + " " + sp[p][0] + " " + sp[p][1] + " " + sp[p][2] + " : ", end='')
				uup = "%s %s" % (is_unique,up)
	#            if is_unique == up == 'U': # uniquely mapped if both ends are
	#                uup = 'U'
	#            else:
	#                uup = 'R'
				print("# " + str(mapq) + " " + str(mp) + " : " + uup) # old: str((mapq+mp)/2.0) + " " + uup
				# note: replace avg of mapq's and display both

			# get new end for next iter
			e = f.readline()
			if not e:
				break

			name, count, mapq, is_unique, snps = parse_line(e)
		else:
			if count > 1: # simply print end
				for p in sorted(snps.keys()) :
					print(p + " " + snps[p][0] + " " + snps[p][1] + " " + snps[p][2] + " : ", end='')
				print("# " + str(mapq) + " : " + is_unique)
			else:
				print('not printing', snps, file=sys.stderr)
			e = ep # and use ep for end of next iter
			name = np
			count = cp
			mapq = mp
			snps = sp
			is_unique = up


if __name__ == '__main__':
	main()
