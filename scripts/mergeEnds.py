#!/usr/bin/env python

# quick script to merge the (sorted) ends output to get a reads output
# -- Murray Patterson

import sys

def pnc(s) : # print s with no carriage return (needs sys)
    sys.stdout.write(str(s))

def main():
	if len(sys.argv) < 2 :
		print "usage : " + str(sys.argv[0]) + " endsFile"
		sys.exit(0)

	f = open(sys.argv[1],"r")

	# get first end
	e = f.readline()
	if not e :
		print "file is empty"
		sys.exit(0)

	t = e.split()
	n = t[0] # name
	c = 0
	m = 0
	u = 0
	s = {} # snps array
	for i in range(len(t)) :
		if t[i] == ":" :
			s[t[i+1]] = [t[i+2],t[i+3],t[i+4]]
			# store snp position and its [base, allele, quality]
		if t[i] == "#" :
			c = int(t[i+1]) # count
			m = int(t[i+2]) # mapping quality
			u = t[i+3] # unique flag ('U' is for unique, 'R' is for repetitive, adopted from BWA XT tag)

	while 1 :
		ep = f.readline() # get second end
		if not ep : # no ep : end e is unpaired
			if c>1 : # so simply print end e
				for p in sorted(s.keys()) : # careful: the default is lists in no particular order, but we want snps to be ordered on their fragment
					pnc(p + " " + s[p][0] + " " + s[p][1] + " " + s[p][2] + " : ")
				print "# " + str(m) + " : " + u
			break

		tp = ep.split()
		np = tp[0]
		cp = 0
		mp = 0
		up = 0
		sp = {}
		for i in range(len(tp)) :
			if tp[i] == ":" :
				sp[tp[i+1]] = [tp[i+2],tp[i+3],tp[i+4]] # pos and [base,allele,qual]
			if tp[i] == "#" :
				cp = int(tp[i+1]) # count
				mp = int(tp[i+2]) # mapq
				up = tp[i+3] # unique or not, see u above

		if n == np : # end e pairs up with end ep
			# output merged pair (a read)
			uup = 0
			if (c+cp)>1 :
				for p in sorted(s.keys()) :
					pnc(p + " " + s[p][0] + " " + s[p][1] + " " + s[p][2] + " : ")
				pnc("-- : ") # add a symbol for gap in paired-end reads
				for p in sorted(sp.keys()) :
					pnc(p + " " + sp[p][0] + " " + sp[p][1] + " " + sp[p][2] + " : ")
				uup = "%s %s" % (u,up)
	#            if u == up == 'U': # uniquely mapped if both ends are
	#                uup = 'U'
	#            else:
	#                uup = 'R'
				print "# " + str(m) + " " + str(mp) + " : " + uup # old: str((m+mp)/2.0) + " " + uup
				# note: replace avg of mapq's and display both

			# get new end for next iter
			e = f.readline()
			if not e : break
			t = e.split()
			n = t[0]
			c = 0
			m = 0
			u = 0
			s = {}
			for i in range(len(t)) :
				if t[i] == ":" :
					s[t[i+1]] = [t[i+2],t[i+3],t[i+4]]
				if t[i] == "#" :
					c = int(t[i+1])
					m = int(t[i+2])
					u = t[i+3]
		else :
			if c>1 : # simply print end
				for p in sorted(s.keys()) :
					pnc(p + " " + s[p][0] + " " + s[p][1] + " " + s[p][2] + " : ")
				print "# " + str(m) + " : " + u
			e = ep # and use ep for end of next iter
			n = np
			c = cp
			m = mp
			s = sp
			u = up


if __name__ == '__main__':
	main()
