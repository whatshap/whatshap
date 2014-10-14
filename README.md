WhatsHap
========

Compilation/Building
--------------------

Create virtualenv in `venv` subfolder:

    virtualenv -p /usr/bin/python3 venv

Install the requirements:

    venv/bin/pip install -r requirements.txt

Compile C++ sources:

    mkdir build && cd build && cmake ../src && make

Then unpack raw.tgz (make sure that folder `raw` is created) and run snakemake.


To Do/Ideas
-----------

* Unit tests.
* Determine best possible phasing that could be achieved:
    * look at covered heterozygous SNPs
    * use only read pairs that cover at least two heterozygous SNPs
    * find connected components
* Are non-uniquely mapping reads used? (They probably should not be)
* Evaluation. Is this relevant?: http://nar.oxfordjournals.org/content/40/5/2041
* Use record.start (0-based) instead of record.POS
* lines like this are output by phasedblocks.py (note start > stop coord):
scaffold250     Phasing exon    8223    8222    .       +       .       gene_id "8222"; transcript_id "8222.1";
* Remove reads that represent a subset of another one
* Extend the example data set to at least two scaffolds
* Current C++ code uses an iterator over matrix columns as input
* Parallelize by working on multiple chromosomes in parallel
* Determine connected components first, then run algorithm in parallel on all
  of them (need to find conn. comp. also afterwards)

Notes
-----

* There is a step in which variants are re-discovered in the BAM file. This may
  fail when the variant caller has used some type of re-alignment (as
  freebayes does). Would be better to integrate this into the variant caller or
  to get the information out of it. This applies only to indels, which are not
  supported right now anyway.
* Input format for HapCompass: http://www.brown.edu/Research/Istrail_Lab/resources/hapcompass_manual.html#sec11

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
