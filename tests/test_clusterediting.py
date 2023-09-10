import itertools
import math

from whatshap.polyphase.solver import (
    AlleleMatrix,
    ClusterEditingSolver,
    scoreReadset,
    TriangleSparseMatrix,
)
from whatshap.testhelpers import string_to_readset


def test_similarities1():
    reads = """
    001001
    110101
    """
    readset = string_to_readset(reads)
    am = AlleleMatrix(readset)
    similarities = scoreReadset(am, 4, 2, 0.06)

    assert not math.isnan(similarities.get(0, 1))
    assert similarities.get(0, 1) < -6.0


def test_similarities2():
    reads = """
    00000
    00000
    00000
    00000
    11111
    11111
    10101
    10101
    """
    readset = string_to_readset(reads)
    am = AlleleMatrix(readset)
    similarities = scoreReadset(am, 4, 4, 0.06)

    assert similarities.get(0, 1) > 1.0
    assert similarities.get(0, 1) == similarities.get(0, 2) == similarities.get(0, 3)
    assert similarities.get(0, 1) == similarities.get(1, 2) == similarities.get(1, 3)
    assert similarities.get(0, 4) < -8.0
    assert similarities.get(0, 5) < -8.0
    assert similarities.get(0, 6) < -1.0
    assert similarities.get(0, 7) < -1.0
    assert similarities.get(4, 5) > 1.0
    assert similarities.get(4, 6) < -1.0


def test_similarities3():
    reads = """
    00000
    00000
    00000
    00000
    11111
    11111
    10101
    10101
    """
    readset = string_to_readset(reads)
    am = AlleleMatrix(readset)
    similarities = scoreReadset(am, 4, 2, 0.06)

    assert similarities.get(0, 1) > 5.0
    assert similarities.get(0, 1) == similarities.get(0, 2) == similarities.get(0, 3)
    assert similarities.get(0, 1) == similarities.get(1, 2) == similarities.get(1, 3)
    assert similarities.get(0, 4) < -8.0
    assert similarities.get(0, 5) < -8.0
    assert similarities.get(0, 6) < -1.0
    assert similarities.get(0, 7) < -1.0
    assert similarities.get(4, 5) > 1.0
    assert similarities.get(4, 6) > 1.0


def test_similarities4():
    reads = """
    000
     000
      000
    111
     111
      101
     110
    """
    readset = string_to_readset(reads)
    am = AlleleMatrix(readset)
    similarities = scoreReadset(am, 2, 2, 0.06)

    assert similarities.get(0, 1) > 2.0
    assert similarities.get(0, 2) == 0.0
    assert similarities.get(1, 2) > 2.0
    assert similarities.get(0, 3) < -5.0 < similarities.get(1, 3) < 0.0 == similarities.get(2, 3)
    assert similarities.get(4, 6) > similarities.get(4, 5) > 0.0


def test_similarities5():
    reads = """
    000
     000
      000
    111
     111
      101
     110
    """
    readset = string_to_readset(reads)
    am = AlleleMatrix(readset)
    similarities = scoreReadset(am, 2, 3, 0.06)

    assert similarities.get(0, 1) > 1.0
    assert similarities.get(0, 2) == 0.0
    assert similarities.get(1, 2) > 0.5
    assert similarities.get(0, 3) < -5.0 < similarities.get(1, 3) < 0.0 == similarities.get(2, 3)
    assert 0.0 > similarities.get(4, 6) > similarities.get(4, 5)


def test_clusterediting1():
    reads = """
        110000010111
        1100000101
         1000 01
         00 0 0 010
         1000001 11
          1111101
          0 10010 1
           0000 010
           1110
           0000 011
            000  10
            0001011
            0  10110
            00010111
            000 0000
        """

    # construct a ReadSet
    readset = string_to_readset(reads)

    # compute similarities
    am = AlleleMatrix(readset)
    similarities = scoreReadset(am, 3, 3, 0.06)

    # run cluster editing
    clusterediting = ClusterEditingSolver(similarities, False)
    readpartitioning = clusterediting.run()

    print("computed clusters: ", readpartitioning)

    # make sure each read occurs only once
    read_ids = list(itertools.chain.from_iterable(readpartitioning))
    duplicates = set([r for r in read_ids if read_ids.count(r) > 1])
    assert len(duplicates) == 0

    assert any(all(x in c for x in [0, 1, 2, 4, 9, 11, 13]) for c in readpartitioning)
    assert any(all(x in c for x in [3, 7, 10, 14]) for c in readpartitioning)
    assert any(all(x in c for x in [5, 8]) for c in readpartitioning)


def test_clusterediting2():
    reads = """
        000000 00 0 00000 0000 0
             1111 11111
               000 00000 0000000
               111111111
                 1000000000
                  0 00000
                    11111
                    1 1 1111 1111111111
                    111111111111
        """

    # construct a ReadSet
    readset = string_to_readset(reads)

    # compute similarities
    am = AlleleMatrix(readset)
    similarities = scoreReadset(am, 3, 2, 0.06)

    # run cluster editing
    clusterediting = ClusterEditingSolver(similarities, False)
    readpartitioning = clusterediting.run()

    assert any(all(x in c for x in [0, 2, 4, 5]) for c in readpartitioning)
    assert any(all(x in c for x in [1, 3, 6, 7, 8]) for c in readpartitioning)


def test_clusterediting3():
    reads = """
        000000 00 0 00000 0000 0
             1111 11111
               000 00000 0000000
               111111111
                 1000000000
                  0 00000
                    11111
                    1 1 1111 1111111111
                    111111111111
        """

    # construct a ReadSet
    readset = string_to_readset(reads)

    # compute similarities
    am = AlleleMatrix(readset)
    similarities = scoreReadset(am, 3, 2, 0.06)

    # run cluster editing
    clusterediting = ClusterEditingSolver(similarities, False)
    readpartitioning = clusterediting.run()

    assert any(all(x in c for x in [0, 2, 4, 5]) for c in readpartitioning)
    assert any(all(x in c for x in [1, 3, 6, 7, 8]) for c in readpartitioning)


def test_clusterediting4():
    reads = """
    0010111110111111111001111
    111111111111111111111 111
    011011111011111 111001111
    00101 111011111 1110011 1
     11 11111111 111111111111
    1111111111111111111111 11
    0010111110111111111001111
    111111111111111111111 111
    011011111011111 111001111
    011011111011111 111001111
    """
    # construct a ReadSet
    readset = string_to_readset(reads)

    # compute similarities
    am = AlleleMatrix(readset)
    similarities = scoreReadset(am, 5, 3, 0.06)

    # run cluster editing
    clusterediting = ClusterEditingSolver(similarities, False)
    readpartitioning = clusterediting.run()

    assert any(all(x in c for x in [0, 2, 3, 6, 8, 9]) for c in readpartitioning)
    assert any(all(x in c for x in [1, 4, 5, 7]) for c in readpartitioning)


def test_clusterediting5():
    reads = """
    0010111110111111111001111
    111111111111111111111 111
    011011111011111 111001111
    00101 111011111 1110011 1
     11 11111111 111111111111
    1111111111111111111111 11
    0010111110111111111001111
    111111111111111111111 111
    011011111011111 111001111
    011011111011111 111001111
    """
    # construct a ReadSet
    readset = string_to_readset(reads)

    # compute similarities
    am = AlleleMatrix(readset)
    similarities = scoreReadset(am, 5, 3, 0.01)

    # run cluster editing
    clusterediting = ClusterEditingSolver(similarities, False)
    readpartitioning = clusterediting.run()

    assert any(all(x in c for x in [0, 3, 6]) for c in readpartitioning)
    assert any(all(x in c for x in [1, 4, 5, 7]) for c in readpartitioning)
    assert any(all(x in c for x in [2, 8, 9]) for c in readpartitioning)


def test_infinity_edges1():
    sim = TriangleSparseMatrix()
    sim.set(0, 1, 1.0)
    sim.set(0, 2, 2.0)
    sim.set(1, 2, -float("inf"))

    ce = ClusterEditingSolver(sim, False)
    clustering = ce.run()

    assert [0, 2] in clustering
    assert [1] in clustering


def test_infinity_edges2():
    sim = TriangleSparseMatrix()
    sim.set(0, 1, -1.0)
    sim.set(0, 2, -2.0)
    sim.set(1, 2, float("inf"))

    ce = ClusterEditingSolver(sim, False)
    clustering = ce.run()

    assert [1, 2] in clustering
    assert [0] in clustering
