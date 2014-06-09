# kate: syntax Python;

SAMPLE = 'sill2'
CHROMOSOMES = ['scaffold221']
SUBSETS = ['moleculo', 'mp', 'moleculomp']

rule all:
	input: expand('result/{chrom}-{subset}.txt', chrom=CHROMOSOMES, subset=SUBSETS)

rule merge_moleculomp:
	input: bam1='data/{chrom}-moleculo.bam', bam2='data/{chrom}-mp.bam'
	output: bam='data/{chrom}-moleculomp.bam'
	shell:
		'picard-tools MergeSamFiles I={input.bam1} I={input.bam2} O={output.bam}'

rule index_bam:
	input: bam='{name}.bam'
	output: bai='{name}.bam.bai'
	shell:
		'samtools index {input}'

rule getends:
	input: bam='data/{chrom}-{subset}.bam', bai='data/{chrom}-{subset}.bam.bai', vcf='data/{chrom}.vcf'
	output: prewif='tmp/{chrom}-{subset}.pre-wif'
	shell:
		'venv/bin/python scripts/getEnds-vcf.py {input.bam} {input.vcf} {wildcards.chrom} {SAMPLE} > {output.prewif}'

rule sortprewif:
	input: '{f}.pre-wif'
	output: '{f}.sorted-pre-wif'
	shell:
		'sort -k1,1 --stable {input} > {output}'

rule mergeends:
	input: '{f}.sorted-pre-wif'
	output: '{f}.wif'
	shell:
		'/usr/bin/python scripts/mergeEnds.py {input} | sort -k 1,1 -n > {output}'

rule shuffle:
	input: '{f}.wif'
	output: '{f}.shuffled-wif'
	shell:
		'shuf {input} > {output}'

rule slicer:
	input: '{f}.shuffled-wif'
	output: '{f}.slice.00.wif'
	shell:
		'/usr/bin/python scripts/tobis-slicer.py -H 20 {input} {wildcards.f}.slice'

rule dp:
	input: '{f}.slice.00.wif'
	output: '{f}.super-reads.wif'
	shell:
		'build/dp --all_het {input} > {output}'

rule extract_het_pos:
	input: vcf='data/{chrom}.vcf'
	output: 'tmp/{chrom}.positions'
	shell:
		'/usr/bin/python scripts/extract-het-pos.py {wildcards.chrom} {input.vcf} > {output}'

rule superread_to_haplotype:
	input: wif='tmp/{chrom}-{subset}.slice.00.wif', superwif='tmp/{chrom}-{subset}.super-reads.wif', positions='tmp/{chrom}.positions'
	output: 'result/{chrom}-{subset}.txt'
	shell:
		'/usr/bin/python scripts/superread-to-haplotype.py -O {input.wif} {input.superwif} {input.positions} > {output}'
