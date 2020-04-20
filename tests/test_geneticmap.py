import pytest
from whatshap.pedigree import GeneticMapRecombinationCostComputer, ParseError


def test_read_genetic_map(tmp_path):
    path = tmp_path / "genetic.map"
    path.write_text("ignored header\n" "568527 0 0\n" "723891 2.9813105581 0.417644215424158\n")
    _ = GeneticMapRecombinationCostComputer(str(path))


def test_read_wrong_number_of_fields(tmp_path):
    path = tmp_path / "genetic.map"
    path.write_text(
        "ignored header\n" "55550 0 0\n" "568322 0 0 17\n" "723891 2.9813105581 0.417644215424158\n"
    )
    with pytest.raises(ParseError):
        _ = GeneticMapRecombinationCostComputer(str(path))


def test_invalid_int(tmp_path):
    path = tmp_path / "genetic.map"
    path.write_text("ignored header\n" "55550 0 0\n" "abc 0 0\n")
    with pytest.raises(ParseError):
        _ = GeneticMapRecombinationCostComputer(str(path))


def test_invalid_float(tmp_path):
    path = tmp_path / "genetic.map"
    path.write_text("ignored header\n" "55550 0 abc\n")
    with pytest.raises(ParseError):
        _ = GeneticMapRecombinationCostComputer(str(path))
