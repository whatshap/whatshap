import os
from tempfile import TemporaryDirectory

from whatshap.cli.hapcut2vcf import run_hapcut2vcf


def test_hapcut2vcf():
    with TemporaryDirectory() as tempdir:
        out = os.path.join(tempdir, "hapcut.vcf")
        run_hapcut2vcf(
            hapcut="tests/data/pacbio/hapcut.txt", vcf="tests/data/pacbio/variants.vcf", output=out,
        )
