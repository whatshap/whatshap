==========
User guide
==========

Run WhatsHap like this::

    whatshap phase -o phased.vcf input.vcf input.bam

Phasing information is added to the VCF file in a way that is compatible with
GATK’s ReadBackedPhasing. That is, the HP tag denotes which set of phased
variants a variant belongs to. The VCF file can also be gzip-compressed.


Features and limitations
========================

WhatsHap supports phasing of variants in diploid genomes.

Supported variant types are SNVs (single-nucleotide variants), insertions,
deletions, MNPs (multiple adjacent SNVs) and “complex” variants. Complex
variants are those that do not fall in any of the other categories, but
are not structural variants. An example is the variant TGCA → AAC.
Structural variants are not phased.

If no reference sequence is provided (using ``--reference``), only
SNVs, insertions and deletions can be phased.

All variants in the input VCF that are marked as being heterozygous
(genotype 0/1) and that have appropriate coverage are used as input for the core
phasing algorithm. If the algorithm could determine how the variant should be
phased, that information will be added to the variant in the output VCF.

Variants can be left unphased for two reasons: Either the variant type is
not supported or the phasing algorithm could not make a phasing decision.
In both cases, the information from the input VCF is simply copied to output
VCF unchanged.


Recommended workflow
====================

Best phasing results are obtained if you sequence your sample(s) on both PacBio
and Illumina: Illumina for high-quality variant calls and PacBio for its long
reads.

1. Map your reads to the reference, making sure that you assign each read to a
read group (the "@RG" header line in the BAM file). WhatsHap supports VCF files
with multiple samples and in order to determine which reads belong to which
sample, it uses the 'sample name' (SM) of the read group. If you have a single
sample only and no or incorrect read group headers, you can run WhatsHap with
``--ignore-read-groups`` instead.

2. Call variants in your sample(s) using the most accurate reads you have. These
will typically be Illumina reads, resulting in a a set of variant calls you can
be reasonably confident in. If you do not know which variant caller to use, yet,
we recommend FreeBayes, which is fast, Open Source and easy to use. In any case,
you will need a standard VCF file as input for WhatsHap in the next step.

3. Run WhatsHap with the VCF file of high-confidence variant calls (obtained in
the previous step) and with the *longest* reads you have. These will typically
be PacBio reads. Phasing works best with long reads, but WhatsHap can use any
read that covers at least two heterozygous variant calls, so even paired-end or
mate-pair reads are somewhat helpful. If you have multiple sets of reads, you
can combine them by providing multiple BAM files on the command line.


Representation of phasing information in VCFs
=============================================

WhatsHap supports two ways in which it can store phasing information in a VCF
file: The standards-compliant ``PS`` tag and the ``HP`` tag used by GATK’s
ReadBackedPhasing tool. When you run ``whatshap phase``, you can select which
format is used by setting ``--tag=PS`` or ``--tag=HP``.

We will use a small VCF file as an example in the following. Unphased, it
looks like this::

    ##fileformat=VCFv4.1
    #CHROM  POS  ID  REF  ALT  QUAL   FILTER  INFO FORMAT  sample1  sample2
    chr1    100  .   A    T    50.0   .       .    GT      0/1      0/1
    chr1    150  .   C    G    50.0   .       .    GT      0/1      1/1
    chr1    300  .   G    T    50.0   .       .    GT      0/1      0/1
    chr1    350  .   T    A    50.0   .       .    GT      0/1      0/1
    chr1    500  .   A    G    50.0   .       .    GT      0/1      1/1

Note that sample 1 is heterozygous at all shown loci (expressed with
``0/1`` in the ``GT`` field).


Phasing represented by pipe (``|``) notation
--------------------------------------------

The ``GT`` fields can be phased by ordering the alleles by haplotype and
separating them with a pipe symbol (``|``) instead of a slash (``/``)::

    ##fileformat=VCFv4.1
    #CHROM  POS  ID  REF  ALT  QUAL   FILTER  INFO FORMAT  sample1  sample2
    chr1    100  .   A    T    50.0   .       .    GT      0|1      0/1
    chr1    150  .   C    G    50.0   .       .    GT      1|0      0/1
    chr1    300  .   G    T    50.0   .       .    GT      1|0      0/1
    chr1    350  .   T    A    50.0   .       .    GT      0|1      0/1
    chr1    500  .   A    G    50.0   .       .    GT      0|1      1/1

The alleles on one of the haplotypes of sample1 are: A, G, T, T, A.
On the other haplotype, they are: T, C, G, A, G.

Swapping ones and zeros in the ``GT`` fields would result in a VCF file with
the equivalent information.


Phasing represented by PS ("phase set") tag
-------------------------------------------

The pipe notation has problems when not all variants in the VCF file can be
phased. The `VCF specification <https://github.com/samtools/hts-specs>`_
introduces the ``PS`` tag to solve some of them. The ``PS`` is a
unique identifier for a "phase set", which is a set of variants that were
be phased relative to each other. There are usually multiple phase sets in
the file, and variants that belong to the same phase set do not need to
be consecutive in the file::

    ##fileformat=VCFv4.1
    #CHROM  POS  ID  REF  ALT  QUAL   FILTER  INFO FORMAT     sample1      sample2
    chr1    100  .   A    T    50.0   .       .    GT:PS:PQ   0|1:100:22   0/1:.:.
    chr1    150  .   C    G    50.0   .       .    GT:PS:PQ   1|0:100:18   0/1:.:.
    chr1    300  .   G    T    50.0   .       .    GT:PS:PQ   1|0:300:23   0/1:.:.
    chr1    350  .   T    A    50.0   .       .    GT:PS:PQ   0|1:300:42   0/1:.:.
    chr1    500  .   A    G    50.0   .       .    GT:PS:PQ   0|1:100:12   0/1:.:.

This VCF contains two phase sets named ``100`` and ``300``. The names are
arbitrary, but WhatsHap will choose the position of the leftmost variant
of the phase set as its name. The variants at 100, 150 and 500 are in the same
phase set, while the variants at 300 and 350 are in a different phase set.
Such a configuration is typically seen when paired-end or mate-pair reads are
used for phasing.

In the case of WhatsHap, the phase sets are identical to the connected components
of the variant connectivity graph. Two variants in that graph are connected if a
read exists that covers them.

The above example also shows usage of the ``PQ`` tag for "phasing quality".
WhatsHap currently does not add this tag.


Phasing represented by HP tag
-----------------------------

GATK’s ReadBackedPhasing tool uses a different way to represent phased variants.
It is in principle the same as the combination of pipe notation with the ``PS``
tag, but the ``GT`` field is left unchanged and all information is added to a
separate ``HP`` tag ("haplotype identifier") instead. This file encodes the same
information as the example above::

    ##fileformat=VCFv4.1
    #CHROM  POS  ID  REF  ALT  QUAL   FILTER  INFO FORMAT     sample1         sample2
    chr1    100  .   A    T    50.0   .       .    GT:HP      0/1:100-1,100-2      0/1:.:.
    chr1    150  .   C    G    50.0   .       .    GT:HP:PQ   0/1:100-2,100-1:18   0/1:.:.
    chr1    300  .   G    T    50.0   .       .    GT:HP:PQ   0/1:300-2,300-1:23   0/1:.:.
    chr1    350  .   T    A    50.0   .       .    GT:HP:PQ   0/1:300-1,300-2:42   0/1:.:.
    chr1    500  .   A    G    50.0   .       .    GT:HP:PQ   0/1:100-1,100-2:12   0/1:.:.

A few notes:

* ReadBackedPhasing does not add the ``PQ`` to the first variant in a phase set/haplotype
  group. This probably means that the phasing quality is to be interpreted as relative to
  the previous or first variant in the set.
* ReadBackedPhasing does not phase indels
* Discussions on the GATK forum on this topic:
   - https://gatkforums.broadinstitute.org/discussion/4226
   - https://gatkforums.broadinstitute.org/discussion/4038/


Trusting the variant caller
===========================

WhatsHap will trust the variant caller to have made the right decision of
whether a variant is heterozygous or homozygous. If you use the option
``--distrust-genotypes``, then this assumption is softened: An optimal solution
could involve switching a variant from being heterozygous to homozygous.
Currently, if that option is enabled and such a switch occurs, the variant
will simply appear as being unphased. No change of the genotype in the VCF is
done.

If you use this option, fewer variants will be phased.

Note that switching homozygous variants to heterozygous is never possible since
only heterozygous variants are considered for phasing.


.. _phasing-pedigrees:

Phasing pedigrees
=================

WhatsHap can take advantage of pedigree information to obtain a much
better phasing. To turn on pedigree mode, run WhatsHap like this::

	whatshap phase --ped pedigree.ped -o phased.vcf input.vcf input.bam

where ``pedigree.ped`` is a plink-compatible PED file to describe the
relationships between samples and ``input.vcf`` is a multi-sample VCF
with all individuals that should be phased. The reads for all individuals
can be in one or more BAM files. WhatsHap will match them based on sample
names provided in the read groups (just like for the default single-individual
mode).

Phasing in pedigree mode requires costs for recombination events. Per
default, WhatsHap will assume a constant recombination rate across the
chromosome to be phased. The recombination rate (in cM/Mb) can be
changed by providing option ``--recombrate``. The default value of
1.26 cM/Mb is suitable for human genomes.

In order to use region-specific recombination rates, a genetic map file
can be provided via option ``--genmap``. WhatsHap expects a three-column
text file like this::

	position COMBINED_rate(cM/Mb) Genetic_Map(cM)
	55550 0 0
	568322 0 0
	568527 0 0
	721290 2.685807669 0.410292036939447
	723819 2.8222713027 0.417429561063975
	723891 2.9813105581 0.417644215424158
	...

The first (header) line is ignored and the three columns are expected to
give the pysical position (in bp), the local recombination rate between the
given position and the position given in the previous row (in cM/Mb), and
the cumulative genetic distance from the start of the chromosome (in cM).
The above example was taken from the 1000 Genomes genetic map `provided by
SHAPEIT
<https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gmap>`_.
Since genetic map files provide information for only one chromosome, the
``--genmap`` option has to be combined with ``--chromosome``.
