"""
Convert hapCUT output format to VCF

HapCUT’s output is combined with the original VCF and
then written as phased VCF to standard output.

HapCUT 1 and 2 are supported.

HapCUT’s output file format is explained at
<https://github.com/vibansal/hapcut#format-of-input-and-output-files>

HapCUT2’s output format is documented at
<https://github.com/pjedge/hapcut2#output-format>
"""
import logging
import re
import sys
from collections import namedtuple
import itertools
from contextlib import ExitStack

from whatshap.vcf import PhasedVcfWriter
from whatshap.core import Read
from whatshap import __version__

logger = logging.getLogger(__name__)


def add_arguments(parser):
	add = parser.add_argument
	add('-o', '--output', default=sys.stdout,
		help='Output VCF file. If omitted, use standard output.')
	add('vcf', metavar='VCF', help='VCF file')
	add('hapcut', metavar='HAPCUT-RESULT', help='hapCUT result file')


HapCutVariant = namedtuple('HapCutVariant', ['chromosome', 'position', 'haplotype1', 'haplotype2', 'component_id'])


class ParseError(Exception):
	pass


class HapCutParser:
	"""Parse HapCUT results"""

	# Example for HapCUT
	#
	# BLOCK: offset: 282 len: 9 phased: 7 SPAN: 5896 MECscore 0.00 fragments 7
	# 282	0	1	1	1065296	T	C	1|0	1,0:-0.0,-0.0,-1.3:0.0:0.0
	# 283	0	1	1	1066259	G	C	1|0	2,0:-0.1,-0.1,-2.7:0.0:0.0
	# 284	0	1	1	1066282	A	G	1|0	3,0:-0.1,-0.2,-4.2:0.0:0.0
	# ********
	# BLOCK: offset: 291 len: 3 phased: 3 SPAN: 7788 MECscore 0.00 fragments 3
	# 291	1	0	1	1072498	G	C	1|0	0,1:-1.4,-0.0,-0.0:0.0:0.0
	# 292	1	0	1	1077064	C	A	1|0	0,2:-2.6,-0.1,-0.0:0.0:0.0
	# 293	1	0	1	1080286	G	A	1|0	0,3:-4.1,-0.1,-0.1:0.1:0.0
	#
	# HapCUT2 output format is slightly different:
	#
	# BLOCK: offset: 306 len: 37 phased: 27 SPAN: 38662 fragments 140
	# 306     0       1       1       1065296 T       C       0/1     0       0.000000        -0.000000
	# 307     0       1       1       1066259 G       C       0/1     0       0.000000        0.000000

	block_re = re.compile(
		'BLOCK: offset: (?P<offset>\d+) len: (?P<len>\d+) phased: (?P<phased>\d+) SPAN: (?P<span>\d+) '
		'(MECscore (?P<mecscore>\d+\.\d+) )?fragments (?P<fragments>\d+)')

	def __init__(self, file):
		self._file = file

	def __iter__(self):
		"""Yield (chromosome, blocks) pairs"""
		yield from self._by_chromosome()

	def parse_blocks(self):
		"""
		Yield a list of HapCutVariant objects for each connected component ('block')
		"""
		state = 'BLOCK'  # DFA states are BLOCK and VARIANT; they describe what we expect next
		block = []
		for line in self._file:
			if state == 'BLOCK':
				state = 'VARIANT'
				if not line.startswith('BLOCK:'):
					raise ParseError('Expected a new block (line starting with "BLOCK:")')
				m = self.block_re.match(line)
				if not m:
					raise ParseError('BLOCK line malformed')
			elif state == 'VARIANT':
				if line.startswith('********'):
					# End marker reached, yield the list of variants in this connected component
					if block:
						yield block
					state = 'BLOCK'
					block = []
				else:
					fields = line.strip().split()
					if len(fields) not in (9, 11):
						raise ParseError('Expected nine fields (for hapCUT 1) or eleven fields (for hapCUT 2) in variant line')
					variant_id, haplotype_1, haplotype_2, chromosome, position, reference_allele, alternative_allele, genotype = fields[:8]

					if len(fields) == 9:  # hapCUT 1
						# The last fields are not actually used, we just check
						# whether they are formatted correctly
						rest = fields[8]
						fields = rest.split(':')
						if len(fields) == 5:
							if not fields[-1] == 'FV':
								raise ParseError('Expected "FV" after last colon')
							fields = fields[:-1]
						if not len(fields) == 4:
							raise ParseError('Too few elements in last (colon-separated) field')
						# allele_counts, genotype_likelihoods, delta, mec_variant = fields
						# allele_counts = [ int(s) for s in allele_counts.split(',') ]
						# genotype_likelihoods = [ float(s) for s in genotype_likelihoods.split(',') ]
						# delta = float(delta)
						# mec_variant = float(mec_variant)
					elif len(fields) == 11:
						# pruned, switch_qual, flip_qual = fields[8:]
						pass
					if haplotype_1 == '-' or haplotype_2 == '-':
						# This happens in hapCUT 2 sometimes
						continue
					variant_id = int(variant_id)
					haplotype_1 = int(haplotype_1)
					haplotype_2 = int(haplotype_2)
					position = int(position) - 1
					component_id = block[0].position if block else position
					variant = HapCutVariant(chromosome, position, haplotype_1, haplotype_2, component_id)
					block.append(variant)
		if len(block) > 0:
			yield block

	def _by_chromosome(self):
		for chromosome, block in itertools.groupby(self.parse_blocks(), lambda b: b[0].chromosome):
			yield chromosome, list(block)


def run_hapcut2vcf(hapcut, vcf, output=sys.stdout):
	command_line = '(whatshap {}) {}'.format(__version__ , ' '.join(sys.argv[1:]))
	with ExitStack() as stack:
		if isinstance(output, str):
			output = stack.enter_context(open(output, 'w'))

		writer = PhasedVcfWriter(vcf, command_line, out_file=output)
		if len(writer.samples) > 1:
			# This would be easy to support with a --sample command-line parameter,
			# but hapCUT does not seem to support multi-sample VCFs, so something
			# must be wrong anyway.
			logger.error('There is more than one sample in this VCF')
			sys.exit(1)
		sample = writer.samples[0]

		f = stack.enter_context(open(hapcut))
		parser = HapCutParser(f)
		for chromosome, blocks in parser:
			logger.info('Read %d phased blocks for chromosome %s', len(blocks), chromosome)

			# Build one read for each haplotype and the connected components
			haplotypes = [ Read(str(i)) for i in (1, 2)]
			components = dict()
			for block in blocks:
				for variant in block:
					haplotypes[0].add_variant(variant.position, variant.haplotype1, 0)
					haplotypes[1].add_variant(variant.position, variant.haplotype2, 0)
					components[variant.position] = variant.component_id

			sample_superreads = { sample: haplotypes }
			sample_components = { sample: components }
			writer.write(chromosome, sample_superreads, sample_components)


def main(args):
	run_hapcut2vcf(**vars(args))
