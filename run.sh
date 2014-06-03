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
getEnds-vcf.py ${bam} ${vcf} ${chromosome} ${individual} > ${prefix}.pre-wif

# Sort pre-wif file to "pair up" reads
sort -k1,1 -s ${prefix}.pre-wif > ${prefix}.sorted-pre-wif

# Merge lines for two reads in a pair and sort by position: create a WIF
mergeEnds.py ${prefix}.sorted-pre-wif | sort -k 1,1 -g > ${prefix}.wif

# Randomly shuffle WIF file (expected by slicer)
shuf ${prefix}.wif > ${prefix}.shuffled-wif

# Create slices (i.e. "bands" of given max. coverage)
tobis-slicer.py -H ${slice_coverage} ${prefix}.shuffled-wif ${prefix}.slice

# Run WhatsHap core algorithm. Using "--all_het" tells the algorithm
# to "trust" the SNP caller; i.e. assume all considered SNP positions
# to be heterozygous.
dp --all_het ${prefix}.slice.00.wif > ${prefix}.super-reads.wif

# Extract positions of heterozygous SNPS 
# TODO: also treats "1/." as heterozygous right now
extract-het-pos.py ${chromosome} ${vcf} > ${prefix}.positions

# Turn WIF format into a "haplotype string", where 
#  0: ref allele
#  1: alt allele
#  -: unphasable: no coverage of read that covers at least 2 SNPs
#  X: unphasable: there is coverage, but still not phasable (tie)
superread-to-haplotype.py -O ${prefix}.wif ${prefix}.super-reads.wif ${prefix}.positions > ${prefix}.haplotypes
