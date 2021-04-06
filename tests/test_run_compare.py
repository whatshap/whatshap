"""
Tests for 'whatshap compare'
"""

from collections import namedtuple
from whatshap.cli.compare import run_compare, compute_switch_flips_poly, compare_block


def test_compare1(tmp_path):
    outtsv = tmp_path / "output.tsv"
    run_compare(
        vcf=["tests/data/phased1.vcf", "tests/data/phased2.vcf"],
        ploidy=2,
        names="p1,p2",
        tsv_pairwise=outtsv,
        sample="sample1",
    )
    lines = [l.split("\t") for l in open(outtsv)]
    assert len(lines) == 3
    Fields = namedtuple("Fields", [f.strip("#\n") for f in lines[0]])
    entry_chrA, entry_chrB = [Fields(*l) for l in lines[1:]]

    assert entry_chrA.chromosome == "chrA"
    assert entry_chrA.all_assessed_pairs == "4"
    assert entry_chrA.all_switches == "1"
    assert entry_chrA.all_switchflips == "1/0"
    assert entry_chrA.blockwise_hamming == "1"
    assert entry_chrA.largestblock_assessed_pairs == "2"
    assert entry_chrA.largestblock_switches == "1"
    assert entry_chrA.largestblock_hamming == "1"

    assert entry_chrB.chromosome == "chrB"
    assert entry_chrB.all_assessed_pairs == "1"
    assert entry_chrB.all_switches == "0"
    assert entry_chrB.all_switchflips == "0/0"
    assert entry_chrB.blockwise_hamming == "0"
    assert entry_chrB.largestblock_assessed_pairs == "1"
    assert entry_chrB.largestblock_switches == "0"
    assert entry_chrB.largestblock_hamming == "0"


def test_compare2(tmp_path):
    outtsv = tmp_path / "output.tsv"
    run_compare(
        vcf=["tests/data/phased1.vcf", "tests/data/phased2.vcf"],
        ploidy=2,
        names="p1,p2",
        tsv_pairwise=outtsv,
        sample="sample2",
    )
    lines = [l.split("\t") for l in open(outtsv)]
    assert len(lines) == 3
    Fields = namedtuple("Fields", [f.strip("#\n") for f in lines[0]])
    entry_chrA, entry_chrB = [Fields(*l) for l in lines[1:]]

    assert entry_chrA.chromosome == "chrA"
    assert entry_chrA.all_assessed_pairs == "6"
    assert entry_chrA.all_switches == "2"
    assert entry_chrA.all_switchflips == "0/1"
    assert entry_chrA.blockwise_hamming == "1"
    assert entry_chrA.largestblock_assessed_pairs == "5"
    assert entry_chrA.largestblock_switches == "2"
    assert entry_chrA.largestblock_hamming == "1"

    assert entry_chrB.chromosome == "chrB"
    assert entry_chrB.all_assessed_pairs == "1"
    assert entry_chrB.all_switches == "1"
    assert entry_chrB.all_switchflips == "1/0"
    assert entry_chrB.blockwise_hamming == "1"
    assert entry_chrB.largestblock_assessed_pairs == "1"
    assert entry_chrB.largestblock_switches == "1"
    assert entry_chrB.largestblock_hamming == "1"


def test_compare_polyploid1(tmp_path):
    outtsv = tmp_path / "output.tsv"
    run_compare(
        vcf=["tests/data/phased.poly1.vcf", "tests/data/phased.poly2.vcf"],
        ploidy=4,
        names="p1,p2",
        tsv_pairwise=outtsv,
        sample="sample1",
    )
    with open(outtsv) as t:
        lines = [line.split("\t") for line in t]
    assert len(lines) == 3
    Fields = namedtuple("Fields", [f.strip("#\n") for f in lines[0]])
    entry_chr21, entry_chr22 = [Fields(*l) for l in lines[1:]]

    assert entry_chr21.chromosome == "chr21"
    assert entry_chr21.all_assessed_pairs == "1"
    assert entry_chr21.all_switches == "0.0"
    assert entry_chr21.all_switchflips == "0.0/0.0"
    assert entry_chr21.blockwise_hamming == "0.0"
    assert entry_chr21.blockwise_diff_genotypes == "0"
    assert entry_chr21.largestblock_assessed_pairs == "1"
    assert entry_chr21.largestblock_switches == "0.0"
    assert entry_chr21.largestblock_hamming == "0.0"
    assert entry_chr21.largestblock_diff_genotypes == "0"

    assert entry_chr22.chromosome == "chr22"
    assert entry_chr22.all_assessed_pairs == "6"
    assert entry_chr22.all_switches == "1.0"
    assert entry_chr22.all_switchflips == "0.0/0.5"
    assert entry_chr22.blockwise_hamming == "0.5"
    assert entry_chr22.blockwise_diff_genotypes == "0"
    assert entry_chr22.largestblock_assessed_pairs == "5"
    assert entry_chr22.largestblock_switches == "1.0"
    assert entry_chr22.largestblock_hamming == "0.5"
    assert entry_chr22.largestblock_diff_genotypes == "0"


def test_compare_polyploid2(tmp_path):
    outtsv = tmp_path / "output.tsv"
    run_compare(
        vcf=["tests/data/phased.poly1.vcf", "tests/data/phased.poly2.vcf"],
        ploidy=4,
        names="p1,p2",
        tsv_pairwise=outtsv,
        sample="sample2",
    )
    with open(outtsv) as t:
        lines = [line.split("\t") for line in t]
    assert len(lines) == 3
    Fields = namedtuple("Fields", [f.strip("#\n") for f in lines[0]])
    entry_chr21, entry_chr22 = [Fields(*l) for l in lines[1:]]

    assert entry_chr21.chromosome == "chr21"
    assert entry_chr21.all_assessed_pairs == "3"
    assert entry_chr21.all_switches == "0.5"
    assert entry_chr21.all_switchflips in ["0.5/0.0", "0.0/0.5"]
    assert entry_chr21.blockwise_hamming == "0.5"
    assert entry_chr21.blockwise_diff_genotypes == "0"
    assert entry_chr21.largestblock_assessed_pairs == "3"
    assert entry_chr21.largestblock_switches == "0.5"
    assert entry_chr21.largestblock_switchflips in ["0.5/0.0", "0.0/0.5"]
    assert entry_chr21.largestblock_hamming == "0.5"
    assert entry_chr21.largestblock_diff_genotypes == "0"

    assert entry_chr22.chromosome == "chr22"
    assert entry_chr22.all_assessed_pairs == "5"
    assert entry_chr22.all_switches == "1.0"
    assert entry_chr22.all_switchflips in ["1.0/0.0", "0.5/0.5", "0.0/1.0"]
    assert entry_chr22.blockwise_hamming == "1.0"
    assert entry_chr22.blockwise_diff_genotypes == "0"
    assert entry_chr22.largestblock_assessed_pairs == "3"
    assert entry_chr22.largestblock_switches == "0.5"
    assert entry_chr22.largestblock_switchflips in ["0.5/0.0", "0.0/0.5"]
    assert entry_chr22.largestblock_hamming == "0.5"
    assert entry_chr22.largestblock_diff_genotypes == "0"


def test_compare_polyploid3(tmp_path):
    outtsv = tmp_path / "output.tsv"
    run_compare(
        vcf=["tests/data/phased.poly1.vcf", "tests/data/phased.poly3.vcf"],
        ploidy=4,
        names="p1,p2",
        tsv_pairwise=outtsv,
        sample="sample1",
    )
    with open(outtsv) as t:
        lines = [line.split("\t") for line in t]
    assert len(lines) == 3
    Fields = namedtuple("Fields", [f.strip("#\n") for f in lines[0]])
    entry_chr21, entry_chr22 = [Fields(*l) for l in lines[1:]]
    assert entry_chr21.chromosome == "chr21"
    assert entry_chr21.all_assessed_pairs == "2"
    assert entry_chr21.all_switches == "0.0"
    assert entry_chr21.all_switchflips == "0.0/0.0"
    assert entry_chr21.blockwise_hamming == "0.0"
    assert entry_chr21.blockwise_diff_genotypes == "0"
    assert entry_chr21.largestblock_assessed_pairs == "2"
    assert entry_chr21.largestblock_switches == "0.0"
    assert entry_chr21.largestblock_switchflips == "0.0/0.0"
    assert entry_chr21.largestblock_hamming == "0.0"
    assert entry_chr21.largestblock_diff_genotypes == "0"

    assert entry_chr22.chromosome == "chr22"
    assert entry_chr22.all_assessed_pairs == "6"
    assert entry_chr22.all_switches == "0.0"
    assert entry_chr22.all_switchflips == "0.0/0.25"
    assert entry_chr22.blockwise_hamming == "0.25"
    assert entry_chr22.blockwise_diff_genotypes == "1"
    assert entry_chr22.largestblock_assessed_pairs == "4"
    assert entry_chr22.largestblock_switches == "0.0"
    assert entry_chr22.largestblock_switchflips == "0.0/0.25"
    assert entry_chr22.largestblock_hamming == "0.25"
    assert entry_chr22.largestblock_diff_genotypes == "1"


def test_compare_only_snvs(tmp_path):
    outtsv = tmp_path / "output.tsv"
    run_compare(
        vcf=["tests/data/phased1.vcf", "tests/data/phased2.vcf"],
        ploidy=2,
        names="p1,p2",
        tsv_pairwise=outtsv,
        sample="sample2",
        only_snvs=True,
    )
    with open(outtsv) as t:
        lines = [line.split("\t") for line in t]
    assert len(lines) == 3
    Fields = namedtuple("Fields", [f.strip("#\n") for f in lines[0]])
    entry_chrA, entry_chrB = [Fields(*l) for l in lines[1:]]

    assert entry_chrA.chromosome == "chrA"
    assert entry_chrA.all_assessed_pairs == "3"
    assert entry_chrA.all_switches == "2"
    assert entry_chrA.all_switchflips == "0/1"
    assert entry_chrA.largestblock_assessed_pairs == "3"
    assert entry_chrA.largestblock_switches == "2"
    assert entry_chrA.largestblock_hamming == "1"

    assert entry_chrB.chromosome == "chrB"
    assert entry_chrB.all_assessed_pairs == "1"
    assert entry_chrB.all_switches == "1"
    assert entry_chrB.all_switchflips == "1/0"
    assert entry_chrB.largestblock_assessed_pairs == "1"
    assert entry_chrB.largestblock_switches == "1"
    assert entry_chrB.largestblock_hamming == "1"


def test_compare_unphased():
    run_compare(
        vcf=["tests/data/unphased.vcf", "tests/data/unphased.vcf", "tests/data/unphased.vcf"],
        ploidy=2,
        sample="sample1",
    )


def test_compute_switch_flips_poly():
    phasing0 = ["0100", "1011"]
    phasing1 = ["0000", "1111"]
    sfp = compute_switch_flips_poly(phasing0, phasing1, flip_cost=3)
    assert sfp.switches == 2.0
    assert sfp.flips == 0

    phasing = ["00000000", "11111111"]
    truth = ["00000000", "11111111"]
    sfp = compute_switch_flips_poly(phasing, truth)
    assert sfp.flips + sfp.switches == 0.0

    phasing = [[0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1]]
    truth = [[0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 0, 0, 0, 0, 0, 0]]
    sfp = compute_switch_flips_poly(phasing, truth)
    assert sfp.flips + sfp.switches == 0.0

    phasing = [[0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1]]
    truth = [[0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    sfp = compute_switch_flips_poly(phasing, truth)
    assert sfp.flips + sfp.switches == 2.0

    phasing = [[1, 1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1]]
    truth = [[0, 0, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1, 1]]
    sfp = compute_switch_flips_poly(phasing, truth)
    assert sfp.flips + sfp.switches == 1.0

    phasing = [[1, 1, 1, 1, 0, 0, 1, 0], [0, 0, 0, 0, 1, 1, 1, 1]]
    truth = [[0, 0, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1, 1]]
    sfp = compute_switch_flips_poly(phasing, truth)
    assert sfp.flips + sfp.switches == 1.5

    phasing = [[1, 1, 1, 1, 0, 0, 1, 0], [0, 0, 0, 0, 1, 1, 1, 1]]
    truth = [[0, 0, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1, 1]]
    sfp = compute_switch_flips_poly(phasing, truth, flip_cost=5, switch_cost=1)
    assert sfp.flips * 5 + sfp.switches == 3.5

    phasing = [[1, 1, 1, 1, 0, 0, 1, 0], [0, 0, 0, 0, 1, 1, 1, 1]]
    truth = [[0, 0, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1, 1]]
    sfp = compute_switch_flips_poly(phasing, truth, flip_cost=1, switch_cost=10)
    assert sfp.flips + sfp.switches * 10 == 3.5

    phasing = [[0, 0, 0, 1, 0, 0, 0, 0], [1, 1, 1, 0, 1, 1, 1, 1]]
    truth = [[0, 0, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1, 1]]
    sfp = compute_switch_flips_poly(phasing, truth)
    assert sfp.flips + sfp.switches == 1.0

    phasing = [[0, 0, 0, 1, 0, 0, 0, 0], [1, 1, 1, 0, 1, 1, 1, 1]]
    truth = [[0, 0, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1, 1]]
    sfp = compute_switch_flips_poly(phasing, truth, flip_cost=5, switch_cost=1)
    assert sfp.flips * 5 + sfp.switches == 2.0

    phasing = [[0, 0, 0, 1, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1, 1]]
    truth = [[0, 0, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1, 1]]
    sfp = compute_switch_flips_poly(phasing, truth, flip_cost=float("inf"), switch_cost=1)
    assert sfp.flips * float("inf") + sfp.switches == float("inf")


def test_compare_block():
    phasing = ["1111111111", "0000000000"]
    truth = ["1111100000", "0000011111"]
    phasing_errors = compare_block(phasing, truth)
    assert phasing_errors.switches == 1
    assert phasing_errors.hamming == 5

    phasing = ["000000", "101111", "111010"]
    truth = ["000000", "101010", "111111"]
    phasing_errors = compare_block(phasing, truth)
    assert phasing_errors.hamming == 2.0 / 3.0
    switch_flips = phasing_errors.switch_flips
    assert switch_flips.switches == 2.0 / 3.0

    phasing = ["1110001", "1011101", "0000010"]
    truth = ["1110001", "1010010", "0001101"]
    phasing_errors = compare_block(phasing, truth)
    assert phasing_errors.hamming == 4.0 / 3.0
    switch_flips = phasing_errors.switch_flips
    assert switch_flips.switches == 2.0 / 3.0

    phasing = ["1111101", "1010001", "0000010"]
    truth = ["1110001", "1010010", "0001101"]
    phasing_errors = compare_block(phasing, truth)
    assert phasing_errors.hamming == 6.0 / 3.0
    switch_flips = phasing_errors.switch_flips
    assert switch_flips.switches == 3.0 / 3.0

    phasing = ["111111", "111111", "111111"]
    truth = ["111111", "000000", "111111"]
    phasing_errors = compare_block(phasing, truth)
    assert phasing_errors.hamming == 2.0
    switch_flips = phasing_errors.switch_flips
    assert switch_flips.switches == 0.0
