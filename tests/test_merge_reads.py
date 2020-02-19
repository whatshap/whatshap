from whatshap.merge import ReadMerger
from whatshap.testhelpers import string_to_readset


def assert_variants(reads, expected):
    # assert that the lists of variants (pos, allele, weight) are identical
    for read, expected_read in zip(reads, expected):
        assert list(read) == list(expected_read)


def test_read_merging():
    reads = string_to_readset(
        """
      0 000000
      111
      11 00111101
      0 00000
    """,
        """
      1 523428
      714
      86 03158958
      8 46626
    """,
    )

    merger = ReadMerger(0.15, 0.25, 100000, 1000)
    merged_reads = merger.merge(reads)
    # default parameter settings

    expected = string_to_readset(
        """
      0 000000
      111
      11 00111101
    """,
        """
      9 989688
      714
      86 03158958
    """,
    )

    assert_variants(merged_reads, expected)


def test_read_merging2():
    reads = string_to_readset(
        """
      0 000000
      111
      11 00111101
      0 00000
    """,
        """
      1 523428
      714
      86 03158958
      8 46626
    """,
    )
    merger = ReadMerger(0.5, 0.5, 1000, 100000)
    merged_reads = merger.merge(reads)
    # error rates and thresholds so high that no merging occurs

    assert_variants(merged_reads, reads)
