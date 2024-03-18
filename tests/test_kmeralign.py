"""
Tests for kmer-alignment performed by kmerald
"""

from whatshap.align import enumerate_all_kmers
from whatshap.align import kmer_align


non_enumerable = [b"", b"A", b"AC", b"ATCGC", b"NNNNNNNNNNN", b"BANANA", b"MISSISSIPPI"]

costs = {
    (53, 214): 2,
    (53, 858): 1,
    (53, 362): 0.5,
    (215, 53): 0.25,
    (215, 214): 1,
    (215, 858): 0.5,
    (215, 362): 0.25,
    (862, 53): 1,
    (862, -5): 10,
    (862, 362): 2.5,
    (378, 53): 0.1,
    (378, 214): 0.25,
    (378, -5): 5,
}
gap = 5
SEQ_1 = [b"AAACCCG", b"AAACCCGG", b"AAATTTCCCG", b"AAAA", b""]
SEQ_2 = [b"GCCCAAA", b"GGCCCAAA", b"GCCCTTTAAA", b"AAAA", b""]


def test_enumeration():
    for string in non_enumerable:
        assert len(enumerate_all_kmers(string, 9)) == 0
        assert len(enumerate_all_kmers(string, 7)) == 0
        assert len(enumerate_all_kmers(string, 6)) == 0
    assert list(enumerate_all_kmers(b"TAAATCCTGG", 7)) == [12341, 215, 862, 3450]
    assert list(enumerate_all_kmers(b"TAAATCCTGG", 11)) == []


def test_kmeralign():
    seq1 = enumerate_all_kmers(b"AATCCTGG", 5)
    seq2 = enumerate_all_kmers(b"AATCCGGG", 5)
    assert kmer_align(seq1, seq2, costs, 5) == 13
    assert kmer_align(seq2, seq1, costs, 5) == 30
    for s1 in SEQ_1:
        for s2 in SEQ_2:
            e_s1 = enumerate_all_kmers(s1, 5)
            e_s2 = enumerate_all_kmers(s2, 5)
            if e_s1 != e_s2:
                expected = gap * (len(e_s1) + len(e_s2))
            else:
                expected = 0
            assert kmer_align(e_s1, e_s2, costs, gap) == expected
