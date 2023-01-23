from whatshap.polyphase import PolyphaseBlockResult
from whatshap.polyphase.solver import AlleleMatrix
from whatshap.polyphase.reorder import find_subinstances, integrate_sub_results, find_breakpoints
from whatshap.testhelpers import string_to_readset


def get_instance1():
    """
    true haplotypes:
    00101100
    01000101
    10111010
    """
    reads = """
        00101
          10110
            1100
        01000
           0010
            0101
        1011
          1110
            1010
        """
    am = AlleleMatrix(string_to_readset(reads))
    clustering = [[0], [1, 4], [2], [3], [5], [6, 7, 8]]
    threads = [
        [0, 3, 5],
        [0, 3, 5],
        [0, 3, 5],
        [1, 1, 5],
        [1, 1, 5],
        [1, 1, 5],
        [4, 2, 5],
        [4, 2, 5],
    ]
    haplotypes = [[0, 0, 1, 0, 1, 1, 0, 1], [0, 1, 0, 0, 0, 1, 0, 0], [1, 0, 1, 1, 1, 0, 1, 0]]
    return am, clustering, threads, haplotypes


def get_instance2():
    """
    true haplotypes:
    000010000000
    010000020101
    101000101001
    001111111110
    """
    reads = """
    00001
      00100000
          000000
    01000
       000020
          020101
    1010001
       0001
         0101001
    001111
        111111
           11110
    """
    am = AlleleMatrix(string_to_readset(reads))
    clustering = [[0], [1, 4, 7], [2], [3], [5], [6], [8], [9, 10, 11]]
    threads = [
        [0, 3, 5, 7],
        [0, 3, 5, 7],
        [0, 3, 5, 7],
        [1, 1, 1, 7],
        [1, 1, 1, 7],
        [1, 1, 1, 7],
        [1, 6, 1, 7],
        [1, 6, 1, 7],
        [1, 6, 1, 7],
        [4, 6, 2, 7],
        [4, 6, 2, 7],
        [4, 6, 2, 7],
    ]
    haplotypes = [
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1],
        [0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1],
        [1, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0],
        [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
    ]
    return am, clustering, threads, haplotypes


def get_instance3():
    """
    true haplotypes:
    000000000
    110101011
    """
    reads = """
        0000
          0000
           0000
            00000
        1101
          0101
           101
            01011
        """
    am = AlleleMatrix(string_to_readset(reads))
    clustering = [[0, 1, 2, 5, 6], [3], [4], [7]]
    threads = [[0, 2], [0, 2], [0, 0], [0, 0], [0, 0], [0, 0], [0, 3], [1, 3], [1, 3]]
    haplotypes = [[0, 0, 0, 1, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0, 1, 0, 1, 1]]
    return am, clustering, threads, haplotypes


def test_find_breakpoints1():
    am, clustering, threads, haplotypes = get_instance1()
    bp = find_breakpoints(threads)
    assert len(bp) == 2
    assert bp[0].position == 3
    assert bp[1].position == 6
    assert bp[0].haplotypes == bp[1].haplotypes == [0, 1]
    assert bp[0].confidence == bp[1].confidence == 0


def test_find_breakpoints2():
    am, clustering, threads, haplotypes = get_instance2()
    bp = find_breakpoints(threads)
    assert len(bp) == 3
    assert bp[0].position == 3
    assert bp[1].position == 6
    assert bp[2].position == 9
    assert bp[0].haplotypes == bp[1].haplotypes == [0, 1, 2]
    assert bp[2].haplotypes == [0, 2]


def test_find_breakpoints3():
    am, clustering, threads, haplotypes = get_instance3()
    bp = find_breakpoints(threads)
    assert len(bp) == 1
    assert bp[0].position == 6
    assert bp[0].haplotypes == [0, 1]
    assert bp[0].confidence == 0


def test_subinstances1():
    am, clustering, threads, haplotypes = get_instance1()
    sub = find_subinstances(am, clustering, threads, haplotypes)

    assert len(sub) == 1
    assert sub[0][0] == 1
    assert sub[0][1] == [0, 1]
    subm = sub[0][2]
    assert len(subm) == 2
    assert subm.getRead(0) in [[(0, 0)], [(0, 1)]]
    assert subm.getRead(1) in [[(0, 0)], [(0, 1)]]


def test_subinstances2():
    am, clustering, threads, haplotypes = get_instance2()
    sub = find_subinstances(am, clustering, threads, haplotypes)

    assert len(sub) == 2

    assert sub[0][0] == 1
    assert sub[0][1] == [0, 1, 2]
    subm = sub[0][2]
    assert len(subm) == 3
    assert subm.getRead(0) in [[(0, 0)], [(0, 1)]]
    assert subm.getRead(1) in [[(0, 0)], [(0, 1)]]
    assert subm.getRead(2) in [[(0, 0)], [(0, 1)]]

    assert sub[1][0] == 1
    assert sub[1][1] == [0, 2]
    subm = sub[1][2]
    assert len(subm) == 2
    print(subm.getRead(0))
    print(subm.getRead(1))
    assert subm.getRead(0) in [[(0, 0)], [(0, 2)]]
    assert subm.getRead(1) in [[(0, 0)], [(0, 2)]]


def test_subinstances3():
    am, clustering, threads, haplotypes = get_instance3()
    sub = find_subinstances(am, clustering, threads, haplotypes)

    assert len(sub) == 1

    assert sub[0][0] == 0
    assert sub[0][1] == [0, 1]
    subm = sub[0][2]
    assert len(subm) == 5
    assert subm.getRead(0) in [[(0, 0)]]
    for i in range(1, 5):
        assert subm.getRead(i) in [[(0, 0), (1, 0)], [(0, 1), (1, 1)]]


def test_integrate_subresults1():
    am, clustering, threads, haplotypes = get_instance1()
    haplotypes_old = haplotypes[:]
    sub = find_subinstances(am, clustering, threads, haplotypes)
    sub_results = []
    sub_results.append(PolyphaseBlockResult(0, [[0], [1]], [[0, 1]], [[0], [1]], []))
    breakpoints = integrate_sub_results(am, sub, sub_results, threads, haplotypes)
    for bp in breakpoints:
        print(bp.position, bp.haplotypes, bp.confidence)
    assert len(breakpoints) == 2
    assert breakpoints[0].position == 3
    assert breakpoints[1].position == 6
    assert breakpoints[0].haplotypes == breakpoints[1].haplotypes == [0, 1]
    assert haplotypes == haplotypes_old


def test_integrate_subresults2():
    am, clustering, threads, haplotypes = get_instance2()
    haplotypes_old = haplotypes[:]
    sub = find_subinstances(am, clustering, threads, haplotypes)
    sub_results = []
    sub_results.append(PolyphaseBlockResult(0, [[0], [1, 2]], [[0, 1, 1]], [[1], [0], [0]], []))
    sub_results.append(PolyphaseBlockResult(0, [[0], [1]], [[0, 1]], [[0], [2]], []))
    breakpoints = integrate_sub_results(am, sub, sub_results, threads, haplotypes)
    for bp in breakpoints:
        print(bp.position, bp.haplotypes, bp.confidence)
    assert len(breakpoints) == 3
    assert breakpoints[0].position == 3
    assert breakpoints[1].position == 6
    assert breakpoints[2].position == 9
    assert breakpoints[0].haplotypes == breakpoints[1].haplotypes == [0, 1, 2]
    assert breakpoints[2].haplotypes == [0, 2]
    assert haplotypes == haplotypes_old


def test_integrate_subresults3():
    am, clustering, threads, haplotypes = get_instance3()
    sub = find_subinstances(am, clustering, threads, haplotypes)
    sub_results = []
    sub_results.append(
        PolyphaseBlockResult(0, [[0, 1, 2], [3, 4]], [[0, 1], [0, 1]], [[0, 0], [1, 1]], [])
    )
    breakpoints = integrate_sub_results(am, sub, sub_results, threads, haplotypes)
    for bp in breakpoints:
        print(bp.position, bp.haplotypes, bp.confidence)
    assert len(breakpoints) == 1
    assert breakpoints[0].position == 6
    assert breakpoints[0].haplotypes == [0, 1]
    assert haplotypes[0] == [0, 0, 0, 0, 0, 0, 0, 0, 0]
    assert haplotypes[1] == [1, 1, 0, 1, 0, 1, 0, 1, 1]
