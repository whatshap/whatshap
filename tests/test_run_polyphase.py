"""
Integration tests that use the command-line entry point run_polyphase.
"""
import os

from pytest import raises
from whatshap.cli.polyphase import run_polyphase, CommandLineError
from whatshap.vcf import VcfReader


def test_polyphase_short_chr22(tmp_path):
    outvcf = tmp_path / "output.vcf"
    run_polyphase(
        phase_input_files=["tests/data/polyploid.chr22.42M.12k.bam"],
        variant_file="tests/data/polyploid.chr22.42M.12k.vcf",
        ploidy=4,
        ignore_read_groups=True,
        output=outvcf,
    )
    assert os.path.isfile(outvcf)

    tables = list(VcfReader(outvcf, phases=True))
    assert len(tables) == 1
    table = tables[0]
    assert table.chromosome == "chr22"
    assert len(table.variants) == 42
    assert table.samples == ["HG00514_NA19240"]


def test_polyphase_multiple_bam(tmp_path):
    outvcf = tmp_path / "output.vcf"
    run_polyphase(
        phase_input_files=[
            "tests/data/polyploid.human1.chr22.42M.5k.bam",
            "tests/data/polyploid.human2.chr22.42M.5k.bam",
        ],
        variant_file="tests/data/polyploid.multisample.chr22.42M.5k.vcf",
        ploidy=2,
        ignore_read_groups=False,
        output=outvcf,
    )
    assert os.path.isfile(outvcf)

    tables = list(VcfReader(outvcf, phases=True))
    assert len(tables) == 1
    table = tables[0]
    assert table.chromosome == "chr22"
    assert len(table.variants) == 9
    assert set(table.samples) == set(["HG00514", "NA19240"])
    assert not all(p is None for p in table.phases_of("HG00514"))
    assert not all(p is None for p in table.phases_of("NA19240"))


def test_wrong_ploidy(tmp_path):
    outvcf = tmp_path / "output.vcf"
    with raises(CommandLineError):
        run_polyphase(
            phase_input_files=["tests/data/polyploid.chr22.42M.12k.bam"],
            variant_file="tests/data/polyploid.chr22.42M.12k.vcf",
            ploidy=3,
            ignore_read_groups=True,
            output=outvcf,
        )


def test_blockcut_sensitivities(tmp_path):
    """ Ensure that the block cut sets are monotone to the sensitivity"""

    results = []
    for s in range(6):
        outvcf = tmp_path / "output{}.vcf".format(s)
        run_polyphase(
            phase_input_files=["tests/data/polyploid.chr22.42M.12k.bam"],
            variant_file="tests/data/polyploid.chr22.42M.12k.vcf",
            ploidy=4,
            ignore_read_groups=True,
            block_cut_sensitivity=s,
            output=outvcf,
        )
        assert os.path.isfile(outvcf)

        tables = list(VcfReader(outvcf, phases=True))
        assert len(tables) == 1
        block_starts = set(
            [i.block_id for i in tables[0].phases_of("HG00514_NA19240") if i is not None]
        )
        results.append(block_starts)
        print(block_starts)

    for s in range(5):
        assert all(cut in results[s + 1] for cut in results[s])
