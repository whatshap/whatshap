# kate: syntax Python;

GATK = 'java -Xmx6G -jar GenomeAnalysisTK.jar'
PYTHON = 'venv/bin/python'  # Must be Python 3

SAMPLE = 'sill2'
SUBSETS = ['moleculo', 'moleculomp0.1']  # perhaps: moleculomp, mp, mp0.1

"""
Expected input files:
- raw/variants.vcf
- raw/moleculo.bam
- raw/bgi-unfixed.bam
"""


rule all:
	input: expand('result/{subset}.txt', subset=SUBSETS)


rule clean:
	shell:
		"rm -f result/* data/*"

rule symlink:
	input: 'raw/{file}'
	output: 'data/{file,(moleculo.bam|bgi-unfixed.bam|variants.vcf)}'
	shell: 'ln -s ../{input} {output}'

rule rgmp_list:
	'Create list of matepair read groups'
	output: rgs='data/rgs-mp.txt'
	run:
		with open(output.rgs, 'w') as f:
			for rg in '10k 20ka 20kb 20kc 5k bgi2ka bgi2kb'.split():
				print(rg, file=f)

rule subsample_matepairs:
	input: bam='data/matepairs.bam'
	output: bam='data/matepairs{frac}.bam'
	shell:
		"samtools view -b {input.bam} -s {wildcards.frac} > {output.bam}"

rule fix_unmapped_mates:
	"""Mark mates as unmapped if the reference they are mapped to is set to "*".
	"""
	input: bam='data/bgi-unfixed.bam'
	output: bam='data/bgi.bam'
	shell:
		"""samtools view -h {input.bam} | awk -vOFS="\t" '!/^@/ && $7=="*"{{$8=0;$2=or($2,8)}};1' | samtools view -bS - > {output.bam}"""
		#picard-tools FixMateInformation I=scaffold221-fixed-tmp.bam O=scaffold221-fixed.bam

rule mpbam:
	'Create the mate-pair BAM file'
	input: bam='data/bgi.bam', rgs='data/rgs-mp.txt'
	output: bam='data/matepairs.bam'
	run:
		shell('samtools view -b -R {input.rgs} {input.bam} > {output.bam}')

rule merge_moleculomp:
	input: bam1='data/moleculo.bam', bam2='data/matepairs{x}.bam'
	output: bam='data/moleculomp{x}.bam'
	shell:
		'picard-tools MergeSamFiles I={input.bam1} I={input.bam2} O={output.bam}'

rule index_bam:
	input: bam='{name}.bam'
	output: bai='{name}.bam.bai'
	shell:
		'samtools index {input}'

rule whatshap:
	input:
		bam='data/{subset}.bam',
		bai='data/{subset}.bam.bai',
		vcf='data/variants.vcf'
	output:
		txt='result/{subset}.txt'
	shell:
		'{PYTHON} -m whatshap --all-het -H 20 {input.bam} {input.vcf} {SAMPLE} > {output.txt}'


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
		vcf='data/gatkphased-{subset}.vcf',
		idx='data/gatkphased-{subset}.vcf.idx'
	input:
		vcf='data/variants.vcf',
		ref='data/ref.fasta',
		fai='data/ref.fasta.fai',
		dictionary='data/ref.dict',
		bam='data/{subset}.bam',
		bai='data/{subset}.bam.bai'
	log: 'data/gatkphased-{subset}.log'
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
