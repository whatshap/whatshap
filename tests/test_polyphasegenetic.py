from whatshap.cli.polyphasegenetic import determine_pedigree
from whatshap.vcf import VcfReader


def test_read_parent_vcf():
    tables = list(
        VcfReader(
            "tests/data/polyphasegenetic.test.parents.vcf",
            only_snvs=False,
            genotype_likelihoods=False,
            ploidy=4,
            mav=True,
        )
    )
    assert len(tables) == 1
    table = tables[0]

    assert table.chromosome == "StSOLv1.1ch03"
    assert len(table.variants) == 135
    assert table.samples == ["Parent_A", "Parent_B"]

    assert table.variants[0].reference_allele == "C"
    assert table.variants[0].alternative_allele == "T"
    assert table.variants[33].reference_allele == "ACC"
    assert table.variants[33].alternative_alleles == ("A", "AC")
    assert table.variants[34].reference_allele == "C"
    assert table.variants[34].alternative_alleles == ("*", "T")


def test_read_progeny_vcf():
    tables = list(
        VcfReader(
            "tests/data/polyphasegenetic.test.progeny.vcf.gz",
            only_snvs=False,
            genotype_likelihoods=False,
            ploidy=4,
            mav=True,
            allele_depth=True,
        )
    )
    assert len(tables) == 1
    table = tables[0]

    assert table.chromosome == "StSOLv1.1ch03"
    assert len(table.variants) == 198
    assert len(table.samples) == 64
    assert table.samples == ["Progeny_{}".format(i) for i in range(1, 65)]

    assert table.variants[0].reference_allele == "C"
    assert table.variants[0].alternative_allele == "T"
    assert table.variants[19].reference_allele == "T"
    assert table.variants[19].alternative_alleles == ("C", "A")
    assert table.variants[71].reference_allele == "AGT"
    assert table.variants[71].alternative_alleles == ("AGGT", "*")

    ad1 = table.allele_depths_of("Progeny_1")
    ad2 = table.allele_depths_of("Progeny_2")
    ad3 = table.allele_depths_of("Progeny_3")
    ad4 = table.allele_depths_of("Progeny_4")

    assert ad1[0] == (1, 4)
    assert ad2[0] == (9,)
    assert ad3[0] == (7,)
    assert ad4[0] == (4, 3)
    assert ad1[19] == (7,)
    assert table.allele_depths_of("Progeny_62")[44] == (3, 1, 1)


def test_pedigree_1():
    file = "tests/data/polyphasegenetic.ped1.txt"
    samples = ["Parent_A"]
    parents = ["Parent_A", "Parent_B", "p1", "p2", "p3", "p4"]
    samples, sam_to_cop, sam_to_prog = determine_pedigree(file, samples, parents)
    assert sam_to_cop["Parent_A"] == "Parent_B"
    assert sam_to_prog["Parent_A"] == ["p1", "p2", "p3", "p4"]
    assert "Parent_B" not in sam_to_cop
    assert "Parent_B" not in sam_to_prog


def test_pedigree_2():
    file = "tests/data/polyphasegenetic.ped1.txt"
    samples = ["Parent_A", "Parent_B"]
    parents = ["Parent_A", "Parent_B", "p1", "p2", "p3", "p4"]
    samples, sam_to_cop, sam_to_prog = determine_pedigree(file, samples, parents)
    assert sam_to_cop["Parent_B"] == "Parent_A"
    assert sam_to_prog["Parent_B"] == sam_to_prog["Parent_A"] == ["p1", "p2", "p3", "p4"]


def test_pedigree_3():
    file = "tests/data/polyphasegenetic.ped2.txt"
    samples = ["Parent_A"]
    parents = ["Parent_A", "Parent_B", "Parent_C", "Parent_D", "p1", "p2", "p3", "p4", "p5"]
    samples, sam_to_cop, sam_to_prog = determine_pedigree(file, samples, parents)
    assert sam_to_cop["Parent_A"] == "Parent_B"
    assert sam_to_prog["Parent_A"] == ["p1", "p2", "p3"]
    assert "Parent_B" not in sam_to_cop
    assert "Parent_B" not in sam_to_prog
    assert "Parent_C" not in sam_to_cop
    assert "Parent_C" not in sam_to_prog


def test_pedigree_4():
    file = "tests/data/polyphasegenetic.ped2.txt"
    samples = ["Parent_C"]
    parents = ["Parent_A", "Parent_B", "Parent_C", "Parent_D", "p1", "p2", "p3", "p4", "p5"]
    samples, sam_to_cop, sam_to_prog = determine_pedigree(file, samples, parents)
    assert sam_to_cop["Parent_C"] == "Parent_D"
    assert sam_to_prog["Parent_C"] == ["p4", "p5"]


def test_pedigree_5():
    file = "tests/data/polyphasegenetic.ped1.txt"
    samples = ["Parent_A"]
    parents = ["Parent_A", "Parent_B", "p1", "p2"]
    samples, sam_to_cop, sam_to_prog = determine_pedigree(file, samples, parents)
    assert sam_to_cop["Parent_A"] == "Parent_B"
    assert sam_to_prog["Parent_A"] == ["p1", "p2"]
    assert "Parent_B" not in sam_to_cop
    assert "Parent_B" not in sam_to_prog


def test_pedigree_6():
    file = "tests/data/polyphasegenetic.ped1.txt"
    samples = ["Parent_A"]
    parents = ["Parent_A", "Parent_B", "p1", "p2"]
    progeny = ["p3", "p4"]
    samples, sam_to_cop, sam_to_prog = determine_pedigree(file, samples, parents, progeny)
    assert sam_to_cop["Parent_A"] == "Parent_B"
    assert sam_to_prog["Parent_A"] == ["p3", "p4"]
    assert "Parent_B" not in sam_to_cop
    assert "Parent_B" not in sam_to_prog
