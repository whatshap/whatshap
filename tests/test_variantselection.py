from whatshap.cli.polyphasegenetic import PolyphaseGeneticParameter
from whatshap.polyphase.variantselection import VariantInfo, compute_phasable_variants
from whatshap.vcf import VcfReader


def test_variantinfo_1():
    vi = VariantInfo([(1, 0)])
    assert vi.get_phasable() == []

    vi.append("A", "C", 1, 0)
    vi.append("A", "G", 1, 1)
    vi.append("A", "T", 0, 1)
    vi.append("C", "A", 1, 0)
    vi.append("C", "G", 1, 0, skip=True)
    vi.append("C", "C", 1, 0, skip=False)
    assert vi[0].alt_count == 1
    assert vi[1].co_alt_count == 1
    assert vi.get_phasable() == [0, 3, 5]
    assert vi.get_node_positions() == [0, 3, 5]


def test_variantinfo_2():
    vi = VariantInfo([(1, 0)])
    assert vi.get_phasable() == []

    vi.append("A", "C", 1, 0)
    vi.append("A", "G", 1, 1)
    vi.append("C", "A", 1, 0)
    vi.append("C", "G", 1, 0)

    try:
        vi.remove_phasable(1)
        assert False
    except ValueError:
        pass

    vi.remove_phasable(3)
    assert vi.get_phasable() == [0, 2]


def test_variantinfo_3():
    vi = VariantInfo([(1, 0)])
    assert vi.get_phasable() == []

    vi.append("A", "C", 1, 0)
    vi.append("A", "G", 1, 1)
    vi.append("C", "A", 1, 0)

    assert vi.get_phasable() == [0, 2]
    vi.correct_type(2, 0, 0)
    assert vi.get_phasable() == [0]
    assert vi[2].alt_count == 0


def test_variantinfo_4():
    vi = VariantInfo([(1, 0), (2, 0)])
    assert vi.get_phasable() == []

    vi.append("A", "C", 1, 0)
    vi.append("A", "G", 1, 1)
    vi.append("C", "A", 2, 0)
    vi.append("C", "G", 1, 0)

    assert vi.get_node_positions() == [0, 2, 2, 3]
    vi.correct_type(0, 2, 0)
    assert vi.get_node_positions() == [0, 0, 2, 2, 3]
    vi.correct_type(3, 0, 0)
    assert vi.get_node_positions() == [0, 0, 2, 2]


def test_compute_phasable_variants_1():
    tables = list(
        VcfReader(
            "tests/data/polyphasegenetic.test.parents.vcf",
            only_snvs=False,
            genotype_likelihoods=False,
            ploidy=4,
            mav=True,
        )
    )
    table = tables[0]

    param = PolyphaseGeneticParameter(4, 20, 0.06, 0, 0, True, True, False, "")

    vi = compute_phasable_variants(table, "Parent_A", "Parent_B", param)
    non_phasable = [
        3,
        9,
        10,
        11,
        12,
        20,
        33,
        34,
        36,
        38,
        52,
        53,
        55,
        56,
        57,
        58,
        59,
        60,
        61,
        63,
        64,
        65,
        67,
        68,
        90,
        91,
        92,
        95,
        96,
        99,
        100,
        101,
        102,
        103,
        104,
        105,
        106,
        108,
        109,
        133,
        134,
    ]
    phasable = vi.get_phasable()

    assert len(vi) == 135
    assert all([x not in phasable for x in non_phasable])
    assert [x for x in range(135) if x not in non_phasable] == phasable

    vi = compute_phasable_variants(table, "Parent_B", "Parent_A", param)
    phasable = vi.get_phasable()
    assert phasable == []


def test_compute_phasable_variants_2():
    tables = list(
        VcfReader(
            "tests/data/polyphasegenetic.test.parents.vcf",
            only_snvs=False,
            genotype_likelihoods=False,
            ploidy=4,
            mav=True,
        )
    )
    table = tables[0]

    param = PolyphaseGeneticParameter(4, 20, 0.06, 1, 0, True, True, False, "")

    vi = compute_phasable_variants(table, "Parent_A", "Parent_B", param)
    non_phasable = [
        3,
        9,
        11,
        12,
        20,
        33,
        34,
        36,
        38,
        52,
        53,
        55,
        56,
        57,
        58,
        59,
        60,
        61,
        63,
        64,
        65,
        67,
        68,
        90,
        91,
        92,
        95,
        96,
        99,
        100,
        101,
        102,
        103,
        104,
        105,
        106,
        108,
        109,
        133,
        134,
    ]
    phasable = vi.get_phasable()

    assert len(vi) == 135
    assert all([x not in phasable for x in non_phasable])
    assert [x for x in range(135) if x not in non_phasable] == phasable

    vi = compute_phasable_variants(table, "Parent_B", "Parent_A", param)
    phasable = vi.get_phasable()
    assert phasable == [10]


def test_compute_phasable_variants_3():
    tables = list(
        VcfReader(
            "tests/data/polyphasegenetic.test.parents.vcf",
            only_snvs=False,
            genotype_likelihoods=False,
            ploidy=4,
            mav=True,
        )
    )
    table = tables[0]

    param = PolyphaseGeneticParameter(4, 20, 0.06, 2, 0, True, True, False, "")

    vi = compute_phasable_variants(table, "Parent_A", "Parent_B", param)
    non_phasable = [33, 34, 36, 38, 96, 99, 106]
    phasable = vi.get_phasable()
    print(phasable)

    assert len(vi) == 135
    assert all([x not in phasable for x in non_phasable])
    assert [x for x in range(135) if x not in non_phasable] == phasable

    vi = compute_phasable_variants(table, "Parent_B", "Parent_A", param)
    phasable = vi.get_phasable()
    assert phasable == [10]
