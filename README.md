WhatsHap
========

This repository is for Marcel's changes.


Notes
-----

`sort -g` is like `-n`, but can also sort floating point numbers and knows about
E notation. It can introduce rounding errors and should therefore only be used
when necessary, see <http://stackoverflow.com/questions/1255782/>.


GATK's phased VCF
-----------------

	awk 'length($4) == 1 && length($5) == 1 { print($2, $10)}' scaffold221-gatk-phased.vcf|cut -d: -f1


Herring example
---------------

	mkdir build && cd build && cmake ../src && make && cd ..
	scripts/getEnds-vcf.py scaffold221-moleculo.bam scaffold221.vcf scaffold221 sill2 > scaffold221.pre-wif
	sort -k1,1 --stable scaffold221.pre-wif > scaffold221.sorted-pre-wif
	scripts/mergeEnds.py scaffold221.sorted-pre-wif | sort -k 1,1 -n > scaffold221.wif
	shuf scaffold221.wif > scaffold221.shuffled-wif
	scripts/tobis-slicer.py -H 20 scaffold221.shuffled-wif scaffold221.slice
	build/dp --all_het scaffold221.slice.00.wif > scaffold221.super-reads.wif
	scripts/extract-het-pos.py scaffold221 scaffold221.vcf > scaffold221.positions
	scripts/superread-to-haplotype.py -O scaffold221.wif scaffold221.super-reads.wif scaffold221.positions > result.txt


pre-wif
-------

* one line per read
* for each read, list all the SNPs from the VCF file that intersect the region the read maps to
* for each SNP, put in its
    * position (on the reference)
    * base (nucleotide) that is actually in the read
    * whether the read supports the REF allele (0) or the ALT allele (1)
    * base quality

sorted-pre-wif
--------------

* same as above, but sorted by read name
* stable sort is used, that is, identically named reads are sorted in the same way they were
  output by the getEnds script

wif
---

* contains data from above, but paired-end reads are merged into one line
* single-end reads are output mostly unchanged
* read name is discarded
* only reads with a count of at least 2 (after merging) are kept


result.txt
----------

* two lines are output. they are equivalent: swap 0 and 1 in the first and you get the second.
