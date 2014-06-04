#!/usr/bin/env python
from __future__ import print_function

# quick script to merge the (sorted) ends output to get a reads output
# -- Murray Patterson

import sys


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

	t = e.split()
	n = t[0] # name
	count = 0
	mapq = 0
	u = 0
	snps = {}  # snps array
	for i in range(len(t)) :
		if t[i] == ":" :
			snps[t[i+1]] = [t[i+2], t[i+3], t[i+4]]
			# store snp position and its [base, allele, quality]
		if t[i] == "#" :
			count = int(t[i+1]) # count
			mapq = int(t[i+2]) # mapping quality
			u = t[i+3] # unique flag ('U' is for unique, 'R' is for repetitive, adopted from BWA XT tag)

	while True:
		ep = f.readline() # get second end
		if not ep: # no ep : end e is unpaired
			if count > 1: # so simply print end e
				for p in sorted(snps.keys()) : # careful: the default is lists in no particular order, but we want snps to be ordered on their fragment
					print(p + " " + snps[p][0] + " " + snps[p][1] + " " + snps[p][2] + " : ", end='')
				print("# " + str(mapq) + " : " + u)
			break

		tp = ep.split()
		np = tp[0]
		cp = 0
		mp = 0
		up = 0
		sp = {}
		for i in range(len(tp)) :
			if tp[i] == ":" :
				sp[tp[i+1]] = [tp[i+2], tp[i+3], tp[i+4]] # pos and [base,allele,qual]
			if tp[i] == "#" :
				cp = int(tp[i+1]) # count
				mp = int(tp[i+2]) # mapq
				up = tp[i+3] # unique or not, see u above

		if n == np : # end e pairs up with end ep
			# output merged pair (a read)
			uup = 0
			if (count+cp)>1 :
				for p in sorted(snps.keys()) :
					print(p + " " + snps[p][0] + " " + snps[p][1] + " " + snps[p][2] + " : ", end='')
				print("-- : ", end='') # add a symbol for gap in paired-end reads
				for p in sorted(sp.keys()) :
					print(p + " " + sp[p][0] + " " + sp[p][1] + " " + sp[p][2] + " : ", end='')
				uup = "%s %s" % (u,up)
	#            if u == up == 'U': # uniquely mapped if both ends are
	#                uup = 'U'
	#            else:
	#                uup = 'R'
				print("# " + str(mapq) + " " + str(mp) + " : " + uup) # old: str((mapq+mp)/2.0) + " " + uup
				# note: replace avg of mapq's and display both

			# get new end for next iter
			e = f.readline()
			if not e:
				break
			t = e.split()
			n = t[0]
			count = 0
			mapq = 0
			u = 0
			snps = {}
			for i in range(len(t)):
				if t[i] == ":":
					snps[t[i+1]] = [t[i+2], t[i+3], t[i+4]]
				if t[i] == "#":
					count = int(t[i+1])
					mapq = int(t[i+2])
					u = t[i+3]
		else:
			if count > 1: # simply print end
				for p in sorted(snps.keys()) :
					print(p + " " + snps[p][0] + " " + snps[p][1] + " " + snps[p][2] + " : ", end='')
				print("# " + str(mapq) + " : " + u)
			e = ep # and use ep for end of next iter
			n = np
			count = cp
			mapq = mp
			snps = sp
			u = up


if __name__ == '__main__':
	main()
