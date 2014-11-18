WhatsHap
========

Compilation/Building
--------------------

Create virtualenv in `venv` subfolder:

    virtualenv -p /usr/bin/python3 venv

Install the requirements:

    venv/bin/pip install -r requirements.txt

Compile whatshap:

    venv/bin/python3 setup.py build_ext -i

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

GATK VCF phasing syntax
-----------------------

It adds these format tags:

    ##FORMAT=<ID=HP,Number=.,Type=String,Description="Read-backed phasing haplotype identifiers">
    ##FORMAT=<ID=PQ,Number=1,Type=Float,Description="Read-backed phasing quality">

Example (edited excerpt):

	24  G T 4399.41 GT:AO:DP:GQ:HP:PL:QA:QR:RO      0/1:136:181:99:   24-1,24-2   :4413,0,1289:5040:1568:42
	72  T G 4229.54 GT:AO:DP:GQ:HP:PL:PQ:QA:QR:RO   0/1:133:199:99:   24-1,24-2   :4244,0,1991:  35.77  :4783:2280:65
	84  T G 3027.84 GT:AO:DP:GQ:HP:PL:PQ:QA:QR:RO   0/1:93:181 :99:   24-1,24-2   :3042,0,2873:  98.44  :3391:3203:87
	194 G A  259.80 GT:AO:DP:GQ:HP:PL:PQ:QA:QR:RO   0/1:10:49  :99:   24-1,24-2   :274,0,1205 :  31.77  :354:1389:39
	254 T A 1041.12 GT:AO:DP:GQ:HP:PL:PQ:QA:QR:RO   0/1:31:55  :99:   24-2,24-1   :1055,0,838 :  31.60  :1181:940:24
	448 C T  311.52 GT:AO:DP:GQ:HP:PL:PQ:QA:QR:RO   0/1:12:58  :99:   24-1,24-2   :325,0,1501 :  37.13  :419:1725:46
	653 T G  298.88 GT:AO:DP:GQ:HP:PL:PQ:QA:QR:RO   0/1:9:26   :99:   24-2,24-1   :313,0,587  :  36.98  :358:663:17

* PQ tag is not added for first variant.
* Indels are not phased
* Forum links:
  https://gatkforums.broadinstitute.org/discussion/4226/
  https://gatkforums.broadinstitute.org/discussion/4038/


File formats
============


result.txt
----------

* two lines are output. If option "--all_het" is used when calling dp, they are
  equivalent: swap 0 and 1 in the first and you get the second; if not, some
  positions might be identified as homozygous.


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
