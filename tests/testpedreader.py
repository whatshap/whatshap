import io
from nose.tools import raises
from whatshap.core import NumericSampleIds
from whatshap.pedigree import PedReader, Trio, ParseError


class TestPedReader:
	def test_parse(self):
		numeric_sample_ids = NumericSampleIds()
		assert len(numeric_sample_ids) == 0
		trios = list(PedReader('tests/data/pedigree.ped', numeric_sample_ids))
		assert trios[0] == Trio(child='child1', mother='mother', father='father')
		assert trios[1] == Trio(child='child2', mother='mother', father='father')
		assert trios[2] == Trio(child='father', mother=None, father=None)
		assert trios[3] == Trio(child='mother', mother=None, father=None)
		assert trios[4] == Trio(child='orphan', mother=None, father=None)
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

