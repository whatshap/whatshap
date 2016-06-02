Various notes
=============

* There is a step in which variants are re-discovered in the BAM file. This may
  fail when the variant caller has used some type of re-alignment (as
  freebayes does). Would be better to integrate this into the variant caller or
  to get the information out of it. This applies only to indels, which are not
  supported right now anyway.
* Input format for HapCompass: http://www.brown.edu/Research/Istrail_Lab/resources/hapcompass_manual.html#sec11

File formats
------------

Phasing in VCFs
~~~~~~~~~~~~~~~

* originally only via 0|1 and 1|0 etc per entry
* then a 'phase set' (PS) added to INFO field: entries with same PS are in same
  set of phased genotypes

GATK VCF phasing syntax
~~~~~~~~~~~~~~~~~~~~~~~

It adds these format tags::

    ##FORMAT=<ID=HP,Number=.,Type=String,Description="Read-backed phasing haplotype identifiers">
    ##FORMAT=<ID=PQ,Number=1,Type=Float,Description="Read-backed phasing quality">

Example (edited excerpt)::

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


Allele detection with re-alignment
----------------------------------

WhatsHap can detect which allele a read contains at a variant position by
aligining a part of the read to the two possible haplotypes. The haplotype
for which the alignment is better wins.

Allele detection through re-alignment is enabled when the ``--reference``
parameter is used on the command-line.

Re-alignment in this version detects slightly *fewer* alleles than the old
algorithm, but this is typically justified because the old algorithm gave
wrong results. Re-alignment however correctly detects that both haplotypes are
equally good and then refuses to choose.

The alignment algorithm uses edit distance at the moment, which allows us to
detect alleles correctly most of the time, but does not allow us to make use
of base qualities (in fact, the weighted algorithm degenerates into an
unweighted one). To fix this, we need a better alignment algorithm.

Here are some examples for how re-alignment works.

Insertion next to a SNP
~~~~~~~~~~~~~~~~~~~~~~~

Haplotypes::

    ref:   CCTTAGT
    alt:   CCTCAGT

Alignment as reported in BAM file::

    ref:   CCT-TAGT
    query: CCTCAAGT

The second ``T`` is aligned to an ``A``, which is not one of the expected bases.
Thus, no variant would be detected here.

Re-aligning the query to the "alt" haplotype, we get::

    alt:   CCTCA-GT
    query: CCTCAAGT

This alignment has lower cost and we therefore detect that the allele in this
read is probably the alternative one.


Ambiguous
~~~~~~~~~

This was previously detected incorrectly::

    ref:   TGCTTTAAGG
    alt:   TGCTTTCAGG
    query: TGCCTTCAAGG

Two possible alignments are ::

    ref:   TGC-TTTAAGG
    query: TGCCTTCAAGG

and ::

    alt:   TGCTTTCA-GG
    query: TGCCTTCAAGG

Both have cost two and therefore the correct allele is unclear.
