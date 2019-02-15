"""
Tag reads by haplotype

Sequencing reads are read from file ALIGNMENTS (in BAM format) and tagged reads
are written to stdout.
"""
import logging
import sys
import pysam
from hashlib import md5
from collections import defaultdict

from contextlib import ExitStack
from whatshap import __version__
from whatshap.vcf import VcfReader
from whatshap.core import NumericSampleIds
from whatshap.bam import AlignmentFileNotIndexedError, SampleNotFoundError
from whatshap.timer import StageTimer
from whatshap.variants import ReadSetReader, ReadSetError
from whatshap.utils import IndexedFasta, FastaNotIndexedError


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
	arg('--ignore-linked-read', default=False, action='store_true',
		help='Ignore linkage information stored in BX tags of the reads.')
	arg('--linked-read-distance-cutoff', '-d', metavar='LINKEDREADDISTANCE', default=50000, type=int,
		help='Assume reads with identical BX tags belong to different read clouds if their '
			'distance is larger than LINKEDREADDISTANCE (default: %(default)s).')
	arg('--ignore-read-groups', default=False, action='store_true',
		help='Ignore read groups in BAM/CRAM header and assume all reads come '
		'from the same sample.')
	arg('--sample', dest='given_samples', metavar='SAMPLE', default=[], action='append',
		help='Name of a sample to phase. If not given, all samples in the '
		'input VCF are phased. Can be used multiple times.')
	arg('variant_file', metavar='VCF', help='VCF file with phased variants (must be gzip-compressed and indexed)')
	arg('alignment_file', metavar='ALIGNMENTS',
		help='File (BAM/CRAM) with read alignments to be tagged by haplotype')


def validate(args, parser):
	pass


def md5_of(filename):
	return md5(open(filename,'rb').read()).hexdigest()


def run_haplotag(
		variant_file,
		alignment_file,
		output=None,
		reference=None,
		ignore_linked_read=False,
		given_samples=None,
		linked_read_distance_cutoff=50000,
		ignore_read_groups=False
	):

	timers = StageTimer()
	timers.start('overall')

	with ExitStack() as stack:
		numeric_sample_ids = NumericSampleIds()
		try:
			readset_reader = stack.enter_context(ReadSetReader([alignment_file],
				reference=reference, numeric_sample_ids=numeric_sample_ids, mapq_threshold=0))
		except OSError as e:
			logger.error(e)
			sys.exit(1)
		except AlignmentFileNotIndexedError as e:
			logger.error('The file %r is not indexed. Please create the appropriate BAM/CRAM '
				'index with "samtools index"', str(e))
			sys.exit(1)

		if reference:
			try:
				fasta = stack.enter_context(IndexedFasta(reference))
			except OSError as e:
				logger.error('%s', e)
				sys.exit(1)
			except FastaNotIndexedError as e:
				logger.error('An index file (.fai) for the reference %r could not be found. '
					'Please create one with "samtools faidx".', str(e))
				sys.exit(1)
		else:
			fasta = None

		# require input VCF to be compressed
		if not variant_file.endswith('gz'):
			logger.error('The input VCF must be compressed (vcf.gz).')
			sys.exit(1)

		vcf_reader = VcfReader(variant_file, indels=True, phases=True)
		vcf_samples = set(vcf_reader.samples)
		logger.info('Found %d samples in VCF file', len(vcf_samples))

		# determine which samples to consider
		if ignore_read_groups and not given_samples and len(vcf_reader.samples) > 1:
			logger.error('When using --ignore-read-groups on a VCF with '
					'multiple samples, --sample must also be used.')
			sys.exit(1)

		if not given_samples:
			given_samples = vcf_reader.samples

		for sample in given_samples:
			if sample not in vcf_samples:
				logger.error('Sample %r requested on command-line not found in VCF', sample)
				sys.exit(1)

		# keep only requested samples
		vcf_samples = vcf_samples.intersection(given_samples)

		# determine which samples are in BAM file
		bam_reader = pysam.AlignmentFile(alignment_file)
		read_groups = bam_reader.header.get('RG', []) 
		bam_samples = set( (rg['SM'] if 'SM' in rg else None) for rg in read_groups )
		rg_to_sample = { rg['ID']:rg['SM'] for rg in read_groups if ('ID' in rg) and ('SM' in rg) }
		logger.info('Samples in BAM file: %s', ','.join(bam_samples))
		samples = vcf_samples
		if not ignore_read_groups:
			samples = bam_samples.intersection(vcf_samples)
			if len(samples) == 0:
				logger.error('No common samples between VCF and BAM file. Aborting.')
				sys.exit(1)
			elif len(samples) < len(bam_samples):
				logger.warning('Not adding phase information for sample(s) %s to BAM file, since they are not present in the VCF or were not given using --sample.', ','.join(bam_samples.difference(vcf_samples)))
		else:
			if len(samples) == 0:
				logger.error('No samples present in VCF. In case --sample was used, the requested sample(s) are not present in the VCF. Aborting.')
				sys.exit(1)

		# Prepare header
		# TODO: convince pysam to allow @HS header line
		header = bam_reader.header.to_dict()
		command_line = ' '.join(['whatshap'] + sys.argv[1:])
		PG_entry = {'ID':'whatshap', 'PN':'whatshap', 'VN':__version__, 'CL':command_line, 'm5': md5_of(variant_file)}
		if 'PG' in header:
			header['PG'].append(PG_entry)
		else:
			header['PG'] = [PG_entry]
		if not output:
			output = '-'
		if output.endswith('.cram'):
			# Write CRAM
			if not reference:
				logger.error("When writing CRAM, you need to provide a reference FASTA using --reference")
				#sys.exit(1)
			kwargs = dict(mode='wc', reference_filename=reference)
		else:
			# Write BAM
			kwargs = dict(mode='wb')
		bam_writer = pysam.AlignmentFile(output, header=pysam.AlignmentHeader.from_dict(header), **kwargs)

		chromosome_name = None
		chromosome_id = None
		skipped_vcf_chromosomes = set()
		n_alignments = 0
		n_tagged = 0
		n_multiple_phase_sets = 0

		# map BX tag to assigned haplotype
		BX_tag_to_haplotype = defaultdict(list)

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
					BX_tag_to_haplotype = defaultdict(list)
					logger.info('Processing alignments on chromosome %s', chromosome_name)
					# Read information on this chromsome from VCF
					variant_table = None
					try:
						variant_table = vcf_reader._fetch(chromosome_name)
						logger.info('... found %s variants for chromosome %s in VCF', len(variant_table), chromosome_name)
					except OSError as e:
						logger.error(str(e))
						sys.exit(1)
					except ValueError:
						logger.info('No variants given for chromosome {} in the input VCF.'.format(chromosome_name))

					# maps read name to (haplotype, quality, phaseset)
					read_to_haplotype = {}
					# Read all reads for this chromosome once to create one core.ReadSet per sample
					# this allows to assign phase to paired-end reads based on both reads
					if variant_table is not None:
						for sample in samples:
							genotypes = variant_table.genotypes_of(sample)
							phases = variant_table.phases_of(sample)
							variantpos_to_phaseinfo = {
								v.position:(int(phases[i].block_id),phases[i].phase) for i,v in enumerate(variant_table.variants) if phases[i] is not None
							}
							variants = [
								v for v, gt, phase in zip(variant_table.variants, genotypes, phases) if gt == 1 and phase is not None
							]
							bam_sample = None if ignore_read_groups else sample
							read_set = read_reads(readset_reader, chromosome_name, variants, bam_sample, fasta)

							# map tag --> set of reads
							BX_tag_to_readlist = defaultdict(list)
							for read in read_set:
								if read.has_BX_tag():
									BX_tag_to_readlist[read.BX_tag].append(read)
							# all reads processed so far
							processed_reads = set()
							for read in read_set:
								if read.name in processed_reads:
									continue
								# mapping: phaseset --> phred scaled difference between costs of assigning reads to haplotype 0 or 1
								haplotype_costs = defaultdict(int)
								reads_to_consider = set()

								processed_reads.add(read.name)
								reads_to_consider.add(read)

								# reads with same BX tag need to be considered too (unless --ignore-linked-read is set)
								if read.has_BX_tag() and not ignore_linked_read:
									for r in BX_tag_to_readlist[read.BX_tag]:
										if not r.name in processed_reads:
											# only select reads close to current one
											if abs(read.reference_start - r.reference_start) <= linked_read_distance_cutoff:
												reads_to_consider.add(r)
								for r in reads_to_consider:
									processed_reads.add(r.name)
									for v in r:
										assert v.allele in [0,1]
										phaseset, allele = variantpos_to_phaseinfo[v.position]
										if v.allele == allele:
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
										BX_tag_to_haplotype[read.BX_tag].append((read.reference_start, haplotype, phaseset))
										for r in reads_to_consider:
											read_to_haplotype[r.name] = (haplotype, abs(quality), phaseset)
											logger.debug('Assigned read %s to haplotype %d with a quality of %d based on %d covered variants', r.name, haplotype, quality, len(r))

				# Only attempt to assign phase of neither secondary nor supplementary
				if (not alignment.is_secondary) and (alignment.flag & 2048 == 0):
					try:
						haplotype, quality, phaseset = read_to_haplotype[alignment.query_name]
						alignment.set_tag('HP', haplotype + 1)
						alignment.set_tag('PC', quality)
						alignment.set_tag('PS', phaseset)
						n_tagged += 1
					except KeyError:
						# check if reads with same tag have been assigned
						if alignment.has_tag('BX'):
							read_clouds = BX_tag_to_haplotype[alignment.get_tag('BX')]
							for (reference_start, haplotype, phaseset) in read_clouds:
								if abs(reference_start - alignment.reference_start) <= linked_read_distance_cutoff:
									alignment.set_tag('HP', haplotype + 1)
									alignment.set_tag('PS', phaseset)
									n_tagged += 1
									break
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
