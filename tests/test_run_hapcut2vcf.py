from whatshap.cli.hapcut2vcf import run_hapcut2vcf


def test_hapcut2vcf(tmp_path):
    out = tmp_path / "hapcut.vcf"
    run_hapcut2vcf(
        hapcut="tests/data/pacbio/hapcut.txt", vcf="tests/data/pacbio/variants.vcf", output=out
    )
