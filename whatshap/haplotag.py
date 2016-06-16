"""
Tag reads by haplotype

Sequencing reads are read from file ALIGNMENTS (in BAM format) and tagged reads
are written to stdout.
"""
import logging
import sys
import pysam
import pyfaidx
from hashlib import md5
from collections import defaultdict

from contextlib import ExitStack
from whatshap import __version__
from whatshap.vcf import VcfReader
from whatshap.core import NumericSampleIds
from whatshap.bam import BamIndexingError, SampleNotFoundError
from whatshap.timer import StageTimer
from whatshap.variants import ReadSetReader, ReadSetError

logger = logging.getLogger(__name__)


def read_reads(readset_reader, chromosome, variants, sample, fasta):
	"""Return a sorted ReadSet"""
	logger.info('Detecting alleles in reads mapped to chromosome %s for sample %r ...', chromosome, sample)
	reference = fasta[chromosome] if fasta else None
	try:
		readset = readset_reader.read(chromosome, variants, sample, reference)
	except ReadSetError as e:
		logger.error("%s", e)
		sys.exit(1)

	# TODO is this necessary?
	for read in readset:
		read.sort()
	readset.sort()

	return readset


def add_arguments(parser):
	arg = parser.add_argument
	arg('-o', '--output', default=None,
		help='Output file. If omitted, use standard output.')
	arg('--reference', '-r', metavar='FASTA',
		help='Reference file. Provide this to detect alleles through re-alignment. '
			'If no index (.fai) exists, it will be created')
	arg('variant_file', metavar='VCF', help='VCF file with phased variants (can be gzip-compressed)')
	arg('alignment_file', metavar='ALIGNMENTS', help='File (BAM) with read alignments to be tagged by haplotype')


def validate(args, parser):
	pass


def md5_of(filename):
	return md5(open(filename,'rb').read()).hexdigest()


def run_haplotag(variant_file, alignment_file, output=None, reference=None):
	timers = StageTimer()
	timers.start('overall')

	with ExitStack() as stack:
		numeric_sample_ids = NumericSampleIds()
		try:
			readset_reader = stack.enter_context(ReadSetReader([alignment_file], numeric_sample_ids, mapq_threshold=0))
		except (OSError, BamIndexingError) as e:
			logger.error(e)
			sys.exit(1)
		if reference:
			try:
				fasta = stack.enter_context(pyfaidx.Fasta(reference, as_raw=True))
			except OSError as e:
				logger.error('%s', e)
				sys.exit(1)
		else:
			fasta = None

		vcf_reader = VcfReader(variant_file, indels=True, phases=True)
		vcf_samples = set(vcf_reader.samples)
		logger.info('Found %d samples in VCF file', len(vcf_samples))

		bam_reader = pysam.AlignmentFile(alignment_file)
		read_groups = bam_reader.header.get('RG', []) 
		bam_samples = set( (rg['SM'] if 'SM' in rg else None) for rg in read_groups )
		rg_to_sample = { rg['ID']:rg['SM'] for rg in read_groups if ('ID' in rg) and ('SM' in rg) }
		logger.info('Samples in BAM file: %s', ','.join(bam_samples))
		samples = bam_samples.intersection(vcf_samples)
		if len(samples) == 0:
			logger.error('No common samples between VCF and BAM file. Aborting.')
			sys.exit(1)
		elif len(samples) < len(bam_samples):
			logger.warning('Not adding phase information for sample(s) %s to BAM file, since they are not present in the VCF', ','.join(bam_samples.difference(vcf_samples)))

		# Prepare header
		# TODO: convince pysam to allow @HS header line
		header = bam_reader.header
		command_line = ' '.join(['whatshap'] + sys.argv[1:])
		PG_entry = { 'PN':'whatshap', 'VN':__version__, 'CL':command_line, 'm5': md5_of(variant_file)}
		if 'PG' in header:
			header['PG'].append(PG_entry)
		else:
			header['PG'] = [PG_entry]
		if output:
			bam_writer = pysam.AlignmentFile(output, 'wb', header=header)
		else:
			bam_writer = pysam.AlignmentFile('-', 'wb', header=header)

		chromosome_name = None
		chromosome_id = None
		skipped_vcf_chromosomes = set()
		vcf_iter = iter(vcf_reader)
		n_alignments = 0
		n_tagged = 0
		n_multiple_phase_sets = 0
		for alignment in bam_reader:
			n_alignments += 1
			alignment.set_tag('HP', value=None)
			alignment.set_tag('PC', value=None)
			alignment.set_tag('PS', value=None)
			if not alignment.is_unmapped:
				# Has chromosome changed?
				if chromosome_id != alignment.reference_id:
					chromosome_id = alignment.reference_id
					chromosome_name = alignment.reference_name
					logger.info('Processing alignments on chromosome %s', chromosome_name)
					if chromosome_name in skipped_vcf_chromosomes:
						logger.error('Chromosome records in alignment file and VCF are sorted differently.')
						sys.exit(1)
					# Read information on this chromsome from VCF
					while True:
						variant_table = next(vcf_iter, None)
						if variant_table is None:
							break
						if variant_table.chromosome == chromosome_name:
							logger.info('... found %s variants chromosome %s in VCF', len(variant_table), chromosome_name)
							break
						else:
							skipped_vcf_chromosomes.add(variant_table.chromosome)
					# maps read name to (haplotype, quality, phaseset)
					read_to_haplotype = {}
					# Read all reads for this chromosome once to create one core.ReadSet per sample
					# this allows to assign phase to paired-end reads based on both reads
					if variant_table is not None:
						for sample in samples:
							genotypes = variant_table.genotypes_of(sample)
							phases = variant_table.phases_of(sample)
							variantpos_to_phaseset = {
								v.position:int(phases[i].block_id) for i,v in enumerate(variant_table.variants) if phases[i] is not None
							}
							variants = [
								v for v, gt, phase in zip(variant_table.variants, genotypes, phases) if gt == 1 and phase is not None
							]
							read_set = read_reads(readset_reader, chromosome_name, variants, sample, fasta)
							for read in read_set:
								# mapping: phaseset --> phred scaled difference between costs of assigning reads to haplotype 0 or 1
								haplotype_costs = defaultdict(int)
								for v in read:
									assert v.allele in [0,1]
									phaseset = variantpos_to_phaseset[v.position]
									if v.allele == 0:
										haplotype_costs[phaseset] += v.quality
									else:
										haplotype_costs[phaseset] -= v.quality
								l = list(haplotype_costs.items())
								l.sort(key=lambda t:-abs(t[1]))
								#logger.info('Read %s: %s', read.name, str(l))
								if len(l) > 0:
									if len(l) > 1:
										n_multiple_phase_sets += 1
									phaseset, quality = l[0]
									if quality != 0:
										haplotype = 0 if quality > 0 else 1
										read_to_haplotype[read.name] = (haplotype, abs(quality), phaseset)
										#logger.debug('Assigned read %s to haplotype %d with a quality of %d based on %d covered variants', read.name, haplotype, quality, len(read))

				# Only attempt to assign phase of neither secondary nor supplementary
				if (not alignment.is_secondary) and (alignment.flag & 2048 == 0):
					try:
						haplotype, quality, phaseset = read_to_haplotype[alignment.query_name]
						alignment.set_tag('HP', haplotype + 1)
						alignment.set_tag('PC', quality)
						alignment.set_tag('PS', phaseset)
						n_tagged += 1
					except KeyError:
						pass
			bam_writer.write(alignment)
			if n_alignments % 100000 == 0:
				logger.info('Processed %d alignment records.', n_alignments)

	logger.info('\n== SUMMARY ==')
	logger.info('Total alignments processed:              %12d', n_alignments)
	logger.info('Alignments that could be tagged:         %12d', n_tagged)
	logger.info('Alignments spanning multiple phase sets: %12d', n_multiple_phase_sets)
	bam_writer.close()


def main(args):
	run_haplotag(**vars(args))
