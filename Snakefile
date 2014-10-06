# kate: syntax Python;

GATK = 'java -Xmx6G -jar GenomeAnalysisTK.jar'

SAMPLE = 'sill2'
CHROMOSOMES = ['scaffold221']
SUBSETS = ['moleculo', 'moleculomp0.1']  # perhaps: moleculomp, mp, mp0.1

"""
Expected input files:
- raw/scaffold221.vcf
- raw/scaffold221-moleculo.bam
- raw/scaffold221-unfixed.bam
"""


rule all:
	input: expand('result/{chrom}-{subset}.txt', chrom=CHROMOSOMES, subset=SUBSETS)

rule symlink:
	input: 'raw/scaffold221{file}'
	output: 'data/scaffold221{file,(-moleculo.bam|-unfixed.bam|.vcf)}'
	shell: 'ln -s ../{input} {output}'

rule rgmp_list:
	'Create list of matepair read groups'
	output: rgs='tmp/rgs-mp.txt'
	run:
		with open(output.rgs, 'w') as f:
			for rg in '10k 20ka 20kb 20kc 5k bgi2ka bgi2kb'.split():
				print(rg, file=f)

rule subsample_mp:
	input: bam='data/scaffold221-mp.bam'
	output: bam='data/scaffold221-mp{frac}.bam'
	shell:
		"samtools view -b {input.bam} -s {wildcards.frac} > {output.bam}"

rule fix_unmapped_mates:
	"""Mark mates as unmapped if the reference they are mapped to is set to "*".
	"""
	input: bam='data/scaffold221-unfixed.bam'
	output: bam='data/scaffold221.bam'
	shell:
		"""samtools view -h {input.bam} | awk -vOFS="\t" '!/^@/ && $7=="*"{{$8=0;$2=or($2,8)}};1' | samtools view -bS - > {output.bam}"""
		#picard-tools FixMateInformation I=scaffold221-fixed-tmp.bam O=scaffold221-fixed.bam


rule mpbam:
	'Create the mate-pair BAM file'
	input: bam='data/{chrom}.bam', rgs='tmp/rgs-mp.txt'
	output: bam='data/{chrom}-mp.bam'
	run:
		shell('samtools view -b -R {input.rgs} {input.bam} > {output.bam}')

rule merge_moleculomp:
	input: bam1='data/{chrom}-moleculo.bam', bam2='data/{chrom}-mp{x}.bam'
	output: bam='data/{chrom}-moleculomp{x}.bam'
	shell:
		'picard-tools MergeSamFiles I={input.bam1} I={input.bam2} O={output.bam}'

rule index_bam:
	input: bam='{name}.bam'
	output: bai='{name}.bam.bai'
	shell:
		'samtools index {input}'

rule getends:
	input: bam='data/{chrom}-{subset}.bam', bai='data/{chrom}-{subset}.bam.bai', vcf='data/{chrom}.vcf'
	output: wif='tmp/{chrom}-{subset}.wif'
	shell:
		'venv/bin/python scripts/whatshap.py -H 20 {input.bam} {input.vcf} {wildcards.chrom} {SAMPLE} > {output.wif}'

rule dp:
	input: '{f}.wif'
	output: '{f}.super-reads.wif'
	shell:
		'build/dp --all_het {input} > {output}'

rule superread_to_haplotype:
	input: wif='tmp/{chrom}-{subset}.wif', superwif='tmp/{chrom}-{subset}.super-reads.wif', vcf='data/{chrom}.vcf'
	output: 'result/{chrom}-{subset}.txt'
	shell:
		'venv/bin/python scripts/superread-to-haplotype.py -O {input.wif} {input.superwif} {input.vcf} {wildcards.chrom} > {output}'


## GATK

rule CreateSequenceDict:
	output: '{base}.dict'
	input: '{base}.fasta'
	resources: time=5
	shell:
		"picard-tools CreateSequenceDictionary R={input} O={output}"


rule faidx:
	output: '{base}.fasta.fai'
	input: '{base}.fasta'
	shell:
		"samtools faidx {input}"


rule ReadBackedPhasing:
	output:
		vcf='tmp/gatkphased-{chrom}-{subset}.vcf',
		idx='tmp/gatkphased-{chrom}-{subset}.vcf.idx'
	input:
		vcf='data/{chrom}.vcf',
		ref='data/{chrom}.fasta',
		fai='data/{chrom}.fasta.fai',
		dictionary='data/{chrom}.dict',
		bam='data/{chrom}-{subset}.bam',
		bai='data/{chrom}-{subset}.bam.bai'
	log: 'tmp/gatkphased-{chrom}-{subset}.log'
	shell:
		r"""
		{GATK} \
			-T ReadBackedPhasing \
			-R {input.ref} \
			-I {input.bam} \
			-L {input.vcf} \
			--variant {input.vcf} \
			-o {output.vcf}.incomplete.vcf \
			--maxPhaseSites 10 \
			--phaseQualityThresh 20.0 >& {log} && \
		mv {output.vcf}.incomplete.vcf {output.vcf} && \
		mv {output.vcf}.incomplete.vcf.idx {output.idx}
		"""
