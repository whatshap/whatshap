import io
from pytest import raises
from whatshap.pedigree import PedReader, Trio, ParseError


class TestPedReader:
    def test_parse(self):
        trios = list(PedReader("tests/data/pedigree.ped"))
        assert trios[0] == Trio(child="child1", mother="mother", father="father")
        assert trios[1] == Trio(child="child2", mother="mother", father="father")
        assert trios[2] == Trio(child="father", mother=None, father=None)
        assert trios[3] == Trio(child="mother", mother=None, father=None)
        assert trios[4] == Trio(child="orphan", mother=None, father=None)

    def test_parse_error(self):
        f = io.StringIO("buggy file")
        with raises(ParseError):
            list(PedReader(f))

    def test_duplicate_individual(self):
        f = io.StringIO("f1 c m f 0 1\nf1 c m f 0 1")
        with raises(ParseError):
            list(PedReader(f))
