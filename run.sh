#!/bin/bash

vcf="venter.chr22.phased.vcf"
individual="Venter"
chromosome="chr22"
slice_coverage="20"
bam="Venter-chr22.miseq.length250.cov1.single.bwamem.sorted.bam"
prefix="${bam/\.bam/}"

echo ${prefix}

# Extract reads from (sorted and indexed) BAM files 
# TODO: only uses data that is already phased according to VCF, 
#       i.e. ignores genotypes with "/" (rather than "|")
# this will also pair up paired-end reads and sort the file by position
# Run WhatsHap core algorithm. Using "--all-het" tells the algorithm
# to "trust" the SNP caller; i.e. assume all considered SNP positions
# to be heterozygous.
scripts/whatshap.py --all-het -H ${slice_coverage} ${bam} ${vcf} ${chromosome} ${individual} > ${prefix}.super-reads.wif

# Extract positions of heterozygous SNPS 
# TODO: also treats "1/." as heterozygous right now
extract-het-pos.py ${chromosome} ${vcf} > ${prefix}.positions

# Turn WIF format into a "haplotype string", where 
#  0: ref allele
#  1: alt allele
#  -: unphasable: no coverage of read that covers at least 2 SNPs
#  X: unphasable: there is coverage, but still not phasable (tie)
superread-to-haplotype.py -O ${prefix}.wif ${prefix}.super-reads.wif ${prefix}.positions > ${prefix}.haplotypes
