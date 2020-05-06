from whatshap.cli.split import run_split


def test_split_bam_no_sequence(tmp_path):
    """
    Test that a BAM file w/o sequence records
    can be processed - see issue 215
    """

    expected_output = [
        "205\t1\t0\t0\n",
        "716\t1\t0\t0\n",
        "1613\t0\t0\t1\n",
        "2250\t1\t0\t0\n",
        "3551\t1\t0\t0\n",
        "4385\t1\t0\t0\n",
        "6750\t1\t0\t0\n",
        "11263\t1\t0\t0\n",
        "12930\t0\t1\t0\n",
        "23225\t0\t1\t0\n",
    ]
    rlen_hist = tmp_path / "rlenhist.tsv"
    input_bam = "tests/data/reads-no-sequence.bam"
    input_list = "tests/data/reads-no-sequence.haplotags.tsv"
    run_split(
        input_bam,
        input_list,
        output_h1="/dev/null",
        output_h2="/dev/null",
        output_untagged="/dev/null",
        read_lengths_histogram=rlen_hist,
    )

    with open(rlen_hist, "r") as dump:
        produced_output = dump.readlines()[1:]  # skip header line
        for e, p in zip(expected_output, produced_output):
            assert e == p
