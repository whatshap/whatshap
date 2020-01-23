from whatshap.cli.unphase import run_unphase


def test_unphase(tmpdir):
    out = tmpdir.join("out.vcf")
    run_unphase("tests/data/phased-via-mixed-HP-PS.vcf", str(out))
    with open("tests/data/unphased.vcf") as f:
        expected = f.read()
    assert expected == out.read_text(encoding="ascii")


def test_unphase_string_typed_ps(tmpdir):
    # Ensure a VCF with PS tags of type String (although against VCF spec) can be read
    run_unphase("tests/data/string_typed_ps_tag.vcf", str(tmpdir.join("out.vcf")))
