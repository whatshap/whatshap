==========
User guide
==========

Run WhatsHap like this::

	whatshap -o phased.vcf input.vcf input.bam

Phasing information is added to the VCF file in a way that is compatible with
GATK’s ReadBackedPhasing. That is, the HP tag denotes which set of phased
variants a variant belongs to. The VCF file can also be gzip-compressed.


Features and limitations
========================

WhatsHap currently supports phasing of SNPs (single-nucleotide changes) in
diploid genomes. All SNPs in the input VCF that are marked as being heterozygous
(genotype 0/1) and that have appropriate coverage are used as input for the core
phasing algorithm. If the algorithm could determine how the SNP should be
phased, that information will be added to the variant in the output VCF.

Insertions, deletions and multi-nucleotide changes (MNPs) are not phased, nor
are structural variations. All records that contain such variants will not
be modified by WhatsHap, that is, the records will simply be copied from the
input to the output file.


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

* Input is a multi-sample VCF with all individuals that should be phased
* Use the ``--ped`` option with a plink-compatible PED file to describe the
  relationships between samples.

*
