"""
Remove phasing information from a VCF file

This script removes all types of phasing information from the input VCF and
prints out the modified VCF to standard output. The modifications are:

- The HP, PS and PQ tags are removed
- Phasing in the GT tag (using pipe notation) is removed. The genotypes are
  sorted in ascending order. For example, a GT value of '1|0' is converted
  to '0/1'.

It is not an error if no phasing information was found.
"""
import sys
import logging
import vcf

logger = logging.getLogger(__name__)

TAGS_TO_REMOVE = frozenset(('HP', 'PQ', 'PS'))


def add_arguments(parser):
	add = parser.add_argument
	add('vcf', metavar='VCF', help='VCF file. Use "-" to read from standard input')


def run_unphase(vcf_path, outfile):
	"""
	Read a VCF file, remove phasing information, and write the result to
	outfile, which must be a file-like object.
	"""
	if vcf_path == '-':
		reader = vcf.Reader(fsock=sys.stdin)
	else:
		reader = vcf.Reader(filename=vcf_path)
	if 'phasing' in reader.metadata:
		reader.metadata['phasing'] = []
	for tag in TAGS_TO_REMOVE:
		if tag in reader.formats:
			del reader.formats[tag]
	writer = vcf.Writer(outfile, template=reader)
	for record in reader:
		formats = record.FORMAT.split(':')
		tag_removed = False
		for tag in TAGS_TO_REMOVE:
			if tag in formats:
				formats.remove(tag)
				tag_removed = True
		if tag_removed:
			record.FORMAT = ':'.join(formats)
			if record.FORMAT not in reader._format_cache:
				reader._format_cache[record.FORMAT] = reader._parse_sample_format(record.FORMAT)
		samp_fmt = reader._format_cache[record.FORMAT]

		for call in record.samples:
			if tag_removed or '|' in call.data.GT:
				values = call.data._asdict()
				for tag in TAGS_TO_REMOVE:
					if tag in values:
						del values[tag]
				gt = values['GT']
				if '|' in gt:
					gt_fields = gt.split('|')
					values['GT'] = '/'.join(sorted(gt_fields))
				call.data = samp_fmt(**values)

		writer.write_record(record)


def main(args):
	run_unphase(args.vcf, sys.stdout)
