WhatsHap
========

This repository is for Marcel's changes.

To Do/Ideas
-----------

* Determine best possible phasing that could be achieved:
    * look at covered heterozygous SNPs
    * use only read pairs that cover at least two heterozygous SNPs
    * find connected components
* Are non-uniquely mapping reads used? (They probably should not be)
* Evaluation. Is this relevant?: http://nar.oxfordjournals.org/content/40/5/2041
* Use record.start (0-based) instead of record.POS

Notes
-----

* There is a step in which variants are re-discovered in the BAM file. This may
  fail when the variant caller has used some type of re-alignment (as
  freebayes does). Would be better to integrate this into the variant caller or
  to get the information out of it. This applies only to indels, which are not
  supported right now anyway.

Phasing in VCFs
---------------

* originally only via 0|1 and 1|0 etc per entry
* then a 'phase set' (PS) added to INFO field: entries with same PS are in same
set of phased genotypes
* GATK's ReadBackedPhasing in newer versions use the 'HP' tag instead, see
  https://gatkforums.broadinstitute.org/discussion/4226/
  https://gatkforums.broadinstitute.org/discussion/4038/


GATK's phased VCF
-----------------

	awk 'length($4) == 1 && length($5) == 1 { print($2, $10)}' scaffold221-gatk-phased.vcf|cut -d: -f1


Fix scaffold221
---------------

(perhaps unnecessary)

	samtools view -h scaffold221.bam | awk -vOFS="\t" '!/^@/ && $7=="*"{$8=0;$2=or($2,8)};1'|samtools view -bS - > scaffold221-fixed-tmp.bam
	picard-tools FixMateInformation I=scaffold221-fixed-tmp.bam O=scaffold221-fixed.bam


Herring example
---------------

	mkdir build && cd build && cmake ../src && make && cd ..

	scripts/getEnds-vcf.py data/scaffold221-moleculo.bam data/scaffold221.vcf scaffold221 sill2 > scaffold221.pre-wif
	sort -k1,1 --stable data/scaffold221.pre-wif > data/scaffold221.sorted-pre-wif
	scripts/mergeEnds.py data/scaffold221.sorted-pre-wif | sort -k 1,1 -n > data/scaffold221.wif
	shuf scaffold221.wif > data/scaffold221.shuffled-wif
	scripts/tobis-slicer.py -H 20 data/scaffold221.shuffled-wif data/scaffold221.slice
	build/dp --all_het data/scaffold221.slice.00.wif > data/scaffold221.super-reads.wif
	scripts/extract-het-pos.py scaffold221 data/scaffold221.vcf > data/scaffold221.positions
	scripts/superread-to-haplotype.py -O data/scaffold221.wif data/scaffold221.super-reads.wif data/scaffold221.positions > result.txt


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

* two lines are output. If option "--all_het" is used when calling dp, they are equivalent: swap 0 and 1 in the first and you get the second; if not, some positions might be identified as homozygous.


Experimental setup
------------------

* Call variants (with freebayes) on scaffold221 (1.2 Mbp) of the herring assembly, using all available reads:
	* mate-pair reads (insert sizes 2k, 5k, 10k, 20k)
	* 2x100 paired-end reads (insert sizes 180, 500, 800)
	* 2x150 MiSeq PE reads
	* Moleculo
* Use the following subsets of data for phasing:
	* Moleculo only
	* Moleculo and mate pairs
	* Mate pairs only (?)
* Phase with GATK's ReadBackedPhasing
* Phase with whatshap

	freebayes -m 1 -f scaffold221.fasta scaffold221-all.bam > scaffold221.vcf
