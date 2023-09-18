"""
Tests for kmer-alignment performed by kmerald

"""
from whatshap.align import enumerate_all_kmers as enum
from whatshap.align import kmer_align as k_align


non_enumerable = ["", "A", "AC", "ATCGC", "NNNNNNNNNNN", "BANANA", "MISSISSIPPI"]

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
SEQ_1 = ["AAACCCG", "AAACCCGG", "AAATTTCCCG", "AAAA", ""]
SEQ_2 = ["GCCCAAA", "GGCCCAAA", "GCCCTTTAAA", "AAAA", ""]


def test_enumeration():
    for string in non_enumerable:
        assert len(enum(str(string).encode("UTF-8"), 9)) == 0
        assert len(enum(str(string).encode("UTF-8"), 7)) == 0
        assert len(enum(str(string).encode("UTF-8"), 6)) == 0
    assert list(enum(str("TAAATCCTGG").encode("UTF-8"), 7)) == [12341, 215, 862, 3450]
    assert list(enum(str("TAAATCCTGG").encode("UTF-8"), 11)) == []


def test_kmeralign():
    seq1 = enum(str("AATCCTGG").encode("UTF-8"), 5)
    seq2 = enum(str("AATCCGGG").encode("UTF-8"), 5)
    assert k_align(seq1, seq2, costs, 5) == 13
    assert k_align(seq2, seq1, costs, 5) == 30
    for s1 in SEQ_1:
        for s2 in SEQ_2:
            e_s1 = enum(s1.encode("UTF-8"), 5)
            e_s2 = enum(s2.encode("UTF-8"), 5)
            if e_s1 != e_s2:
                assert k_align(e_s1, e_s2, costs, gap) == gap * (len(e_s1) + len(e_s2))
            else:
                assert k_align(e_s1, e_s2, costs, gap) == 0
