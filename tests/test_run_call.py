from io import StringIO
from whatshap.call import run_call 
import vcf
from pytest import raises
from tempfile import TemporaryDirectory

def test_call():
	with TemporaryDirectory() as tempdir:
		output = tempdir + '/output.vcf'
		run_call('tests/data/pacbio/reference.fasta', 'tests/data/pacbio/pacbio.bam', pacbio=True, output=output)
		computed_lines = []
		expected_lines = []
		for line in open(output, 'r'):
			if line.startswith('#'):
				continue
			computed_lines.append(line)
		for line in open('tests/data/expected-calls.vcf'):
			if line.startswith('#'):
				continue
			expected_lines.append(line)
		assert computed_lines == expected_lines
