from whatshap.cli.find_snv_candidates import run_find_snv_candidates


def test_call(tmpdir):
    output = str(tmpdir.join("output.vcf"))
    run_find_snv_candidates(
        "tests/data/pacbio/reference.fasta",
        "tests/data/pacbio/pacbio.bam",
        datatype="pacbio",
        output=output,
    )
    computed_lines = []
    expected_lines = []
    for line in open(output, "r"):
        if line.startswith("#"):
            continue
        computed_lines.append(line)
    for line in open("tests/data/expected-calls.vcf"):
        if line.startswith("#"):
            continue
        expected_lines.append(line)
    assert computed_lines == expected_lines
