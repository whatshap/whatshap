import gzip
import pyfaidx
import os.path


class UnknownFileFormatError(Exception):
	pass


class FastaNotIndexedError(Exception):
	pass


def detect_file_format(path):
	"""
	Detect file format and return 'BAM', 'CRAM', 'VCF' or raise an UnknownFileFormatError.
	'VCF' is returned for both uncompressed and compressed VCFs (.vcf and .vcf.gz).
	"""
	try:
		with open(path, 'rb') as f:
			first_bytes = f.read(16)
			if first_bytes.startswith(b'CRAM'):
				return 'CRAM'
			if first_bytes.startswith(b'##fileformat=VCF'):
				return 'VCF'

		# Even 'uncompressed' BAM files have a BGZF header and are therefore
		# detected correctly
		with gzip.GzipFile(path, 'rb') as f:
			first_bytes = f.read(16)
			if first_bytes.startswith(b'BAM\1'):
				return 'BAM'
			elif first_bytes.startswith(b'##fileformat=VCF'):
				return 'VCF'
	except OSError:
		pass
	raise UnknownFileFormatError()


def IndexedFasta(path):
	if not os.path.exists(path + '.fai'):
		raise FastaNotIndexedError(path)
	return pyfaidx.Fasta(path, as_raw=True, sequence_always_upper=True)
