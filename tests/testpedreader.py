import io
from nose.tools import raises
from whatshap.core import NumericSampleIds
from whatshap.pedigree import PedReader, Individual, ParseError

class TestPedReader:
	def test_parse(self):
		numeric_sample_ids = NumericSampleIds()
		assert len(numeric_sample_ids) == 0
		individuals = list(PedReader('tests/data/pedigree.ped', numeric_sample_ids))
		assert individuals[0] == Individual(id='child1', mother_id='mother', father_id='father')
		assert individuals[1] == Individual(id='child2', mother_id='mother', father_id='father')
		assert individuals[2] == Individual(id='father', mother_id=None, father_id=None)
		assert individuals[3] == Individual(id='mother', mother_id=None, father_id=None)
		assert individuals[4] == Individual(id='orphan', mother_id=None, father_id=None)
		assert len(numeric_sample_ids) == 5

	@raises(ParseError)
	def test_parse_error(self):
		numeric_sample_ids = NumericSampleIds()
		f = io.StringIO("buggy file")
		list(PedReader(f, numeric_sample_ids))

	@raises(ParseError)
	def test_duplicate_individual(self):
		numeric_sample_ids = NumericSampleIds()
		f = io.StringIO("f1 c m f 0 1\nf1 c m f 0 1")
		list(PedReader(f, numeric_sample_ids))

