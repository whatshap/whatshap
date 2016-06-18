# This Makefile documents how the PacBio test data set was generated.
# The data is from the Genome in a Bottle dataset (individual HG004)

SHELL := /bin/bash

all: reference.fasta pacbio.bam pacbio.bam.bai variants.vcf

ucsc.hg19.fasta:
	@echo "Please obtain file ucsc.hg19.fasta from the GATK resource bundle to continue"
	@exit 1

reference.fasta: ucsc.hg19.fasta ucsc.hg19.fasta.fai
	( echo '>ref'; samtools faidx $< chr6:10029001-10056000 | sed 1d | tr a-z A-Z ) > $@

illumina.bam:
	samtools view -b ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/NIST_Illumina_2x250bps/HG004_250bp_All-29933060/HG004-250bp-All_S1.bam chr6:10039500-10049000 > $@

variants.vcf: ucsc.hg19.fasta illumina.bam illumina.bam.bai
	freebayes -f ucsc.hg19.fasta illumina.bam | sed 's|^chr6|ref|' | awk -vFS='\t' -vOFS='\t' '!/^#/ {$$2-=10029000; $$8="."};1' | grep -v '^##INFO=' | sed -r 's|GT(:[A-Z][A-Z])*\t([01/]*):.*|GT\t\2|' > $@

downloaded.bam:
	samtools view -b ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/PacBio_MtSinai_NIST/CSHL_bwamem_bam_GRCh37/BWA-MEM_Chr6_HG004_merged_11_12.sort.bam 6:10039500-10049000 > $@

pacbio.bam: downloaded.bam reference.fasta
	samtools faidx reference.fasta
	( echo -e "@HD\tVN:1.2\tSO:coordinate\n@SQ\tSN:ref\tLN:27000\n@RG\tID:1\tSM:HG004_250bp_All"; samtools view -s 0.3 -x MD downloaded.bam | awk -vFS='\t' -vOFS='\t' '{$$3="ref"; $$4-=10029000; print $$0 "\tRG:Z:1"}' ) | samtools view -b - > $@

%.bam.bai: %.bam
	samtools index $<

%.fasta.fai: %.fasta
	samtools faidx $<

