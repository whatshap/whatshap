from whatshap.core import (
    ReadSet,
    DynamicSparseGraph,
    ClusterEditingSolver,
    scoreReadsetGlobal,
    scoreReadsetLocal,
)
from whatshap.testhelpers import string_to_readset
import itertools
import math


def test_clusterediting1():

    reads = """
        110000010111
        1100000101
         1000 01
         00 0 0 000
         1000001 11
          1111101
          0 10010 1
           0000 010
           1110
           0000 011
            000  00
            0001011
            0  10110
            00010111
            000 0000
        """

    # construct a ReadSet
    readset = string_to_readset(reads)

    # compute similarities
    similarities = scoreReadsetGlobal(readset, 5, 4)

    # create read graph
    n_reads = len(readset)
    graph = DynamicSparseGraph(n_reads)

    # insert edges
    for (read1, read2) in similarities:
        graph.addEdge(read1, read2, similarities.get(read1, read2))

    # run cluster editing
    clusterediting = ClusterEditingSolver(graph, False)
    readpartitioning = clusterediting.run()

    print("computed clusters: ", readpartitioning)

    # make sure each read occurs only once
    read_ids = list(itertools.chain.from_iterable(readpartitioning))
    duplicates = set([r for r in read_ids if read_ids.count(r) > 1])
    print("duplicates:", duplicates)
    assert len(duplicates) == 0


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
    similarities = scoreReadsetGlobal(readset, 5, 2)
    print(similarities)

    # create read graph
    n_reads = len(readset)
    graph = DynamicSparseGraph(n_reads)

    # insert edges
    for (read1, read2) in similarities:
        graph.addEdge(read1, read2, similarities.get(read1, read2))

    # run cluster editing
    clusterediting = ClusterEditingSolver(graph, False)
    readpartitioning = clusterediting.run()

    print("computed clusters: ", readpartitioning)

    # make sure each read occurs only once
    read_ids = list(itertools.chain.from_iterable(readpartitioning))
    duplicates = set([r for r in read_ids if read_ids.count(r) > 1])
    print("duplicates:", duplicates)
    assert len(duplicates) == 0


def test_clusterediting3():
    reads = """
    0010111110111111111001111
    111111111111111111111 111
    011011111011111 111001111
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
    similarities = scoreReadsetGlobal(readset, 5, 3)

    print(similarities)

    # create read graph
    n_reads = len(readset)
    graph = DynamicSparseGraph(n_reads)

    # insert edges
    for (read1, read2) in similarities:
        graph.addEdge(read1, read2, similarities.get(read1, read2))

    # run cluster editing
    clusterediting = ClusterEditingSolver(graph, False)
    readpartitioning = clusterediting.run()

    print("computed clusters: ", readpartitioning)


def test_similarities1():
    reads = """
    001001
    110101
    """
    readset = string_to_readset(reads)
    similarities = scoreReadsetGlobal(readset, 4, 2)
    # computed similarity is 'nan'
    print("computed similarities:", similarities)
    assert not math.isnan(similarities.get(0, 1))


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
    similarities = scoreReadsetGlobal(readset, 4, 4)
    print("computed similarities:", similarities)
