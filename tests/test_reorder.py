from whatshap.polyphase import PolyphaseBlockResult, PhaseBreakpoint
from whatshap.polyphase.solver import AlleleMatrix
from whatshap.polyphase.reorder import (
    find_subinstances,
    integrate_sub_results,
    find_breakpoints,
    get_heterozygous_pos_for_haps,
    compute_link_likelihoods,
    compute_phase_affiliation,
    get_optimal_assignments,
)
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


def get_test_instance4():
    return [
        [0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0],
        [0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0],
    ]


def test_get_heterozygous_pos_for_haps1():
    haplotypes = get_test_instance4()
    l, r = get_heterozygous_pos_for_haps(haplotypes, [0, 1], 6, limit=1)
    assert l == [3]
    assert r == [7]
    l, r = get_heterozygous_pos_for_haps(haplotypes, [0, 1], 6, limit=2)
    assert l == [2, 3]
    assert r == [7, 9]


def test_get_heterozygous_pos_for_haps2():
    haplotypes = get_test_instance4()
    l, r = get_heterozygous_pos_for_haps(haplotypes, [0, 1], 7, limit=2)
    assert l == [2, 3]
    assert r == [7, 9]
    l, r = get_heterozygous_pos_for_haps(haplotypes, [0, 1], 7, limit=3)
    assert l == [2, 3]
    assert r == [7, 9]


def test_get_heterozygous_pos_for_haps3():
    haplotypes = get_test_instance4()
    l, r = get_heterozygous_pos_for_haps(haplotypes, [0, 2], 3, limit=2)
    assert l == [1]
    assert r == []
    l, r = get_heterozygous_pos_for_haps(haplotypes, [0, 1, 2], 3, limit=2)
    assert l == [1, 2]
    assert r == [3, 7]


def test_compute_link_likelihoods():
    am, clustering, threads, haplotypes = get_instance2()
    bp = [
        PhaseBreakpoint(3, [0, 1, 2], 0),
        PhaseBreakpoint(6, [0, 1, 2], 0),
        PhaseBreakpoint(9, [0, 1], 0),
    ]
    llh = compute_link_likelihoods(threads, haplotypes, bp, clustering, am, 0.07)
    assert llh[0][(0, 2, 1)] > llh[0][(0, 1, 2)]
    assert llh[0][(1, 0, 2)] <= llh[0][(0, 1, 2)]
    assert llh[0][(1, 2, 0)] < llh[0][(0, 2, 1)]
    assert llh[0][(2, 0, 1)] < llh[0][(0, 2, 1)]
    assert llh[0][(2, 1, 0)] < llh[0][(0, 1, 2)]

    assert llh[1][(0, 1, 2)] == max(llh[1].values())

    assert llh[2][(0, 1)] == max(llh[2].values())


def test_compute_phase_affiliation():
    am, clustering, threads, haplotypes = get_instance2()
    bp = [
        PhaseBreakpoint(3, [0, 1, 2], 0),
        PhaseBreakpoint(6, [0, 1, 2], 0),
        PhaseBreakpoint(9, [0, 1], 0),
    ]
    superreads = """
    0  01  0   0
    0  00  2   1
    1  00  0   1
    0  11  1   0
    """
    pp = AlleleMatrix(string_to_readset(superreads))
    aff = compute_phase_affiliation(am, haplotypes, bp, pp, 0.07)

    assert len(aff) == 4
    assert aff[0][0][0] == max(aff[0][0])
    assert aff[0][1][1] == max(aff[0][1])
    assert aff[0][2][2] == max(aff[0][2])
    assert aff[0][3][3] == max(aff[0][3])

    assert aff[1][0][0] == max(aff[1][0])
    assert aff[1][1][0] == max(aff[1][1])
    assert aff[1][2][2] == max(aff[1][2])
    assert aff[1][3][3] == max(aff[1][3])

    assert aff[2][0][0] == max(aff[2][0])
    assert aff[2][1][2] == max(aff[2][1])
    assert aff[2][2][1] == max(aff[2][2])
    assert aff[2][3][3] == max(aff[2][3])

    assert aff[3][0][1] == max(aff[3][0])
    assert aff[3][1][2] == max(aff[3][1])
    assert aff[3][2][0] == max(aff[3][2])
    assert aff[3][3][3] == max(aff[3][3])


def test_get_optimal_permutations1():
    am, clustering, threads, haplotypes = get_instance2()
    bp = [
        PhaseBreakpoint(3, [0, 1, 2], 0),
        PhaseBreakpoint(6, [0, 1, 2], 0),
        PhaseBreakpoint(9, [0, 1], 0),
    ]
    lllh = compute_link_likelihoods(threads, haplotypes, bp, clustering, am, 0.07)
    asmnts = get_optimal_assignments(bp, lllh, 4, None)
    assert asmnts[0] == [0, 1, 2, 3]
    assert asmnts[1] in [[0, 1, 2, 3], [0, 2, 1, 3], [1, 0, 2, 3], [2, 0, 1, 3]]
    assert (asmnts[2] in [[0, 2, 1, 3], [2, 0, 1, 3]]) or (
        asmnts[3] in [[1, 2, 0, 3], [1, 2, 3, 0], [2, 1, 0, 3], [2, 1, 3, 0]]
    )
    assert asmnts[2][2:] == asmnts[3][2:]


def test_get_optimal_permutations2():
    am, clustering, threads, haplotypes = get_instance2()
    bp = [
        PhaseBreakpoint(3, [0, 1, 2], 0),
        PhaseBreakpoint(6, [0, 1, 2], 0),
        PhaseBreakpoint(9, [0, 1], 0),
    ]
    lllh = compute_link_likelihoods(threads, haplotypes, bp, clustering, am, 0.07)
    superreads = """
    0  01  0   0
    0  00  2   1
    1  00  0   1
    0  11  1   0
    """

    pp = AlleleMatrix(string_to_readset(superreads))
    aff = compute_phase_affiliation(am, haplotypes, bp, pp, 0.07)
    asmnts = get_optimal_assignments(bp, lllh, 4, aff)
    assert asmnts[0] == [0, 1, 2, 3]
    assert asmnts[1] in [[0, 1, 2, 3], [0, 2, 1, 3], [1, 0, 2, 3], [2, 0, 1, 3]]
    assert (asmnts[2] in [[0, 2, 1, 3], [2, 0, 1, 3]]) or (
        asmnts[3] in [[1, 2, 0, 3], [1, 2, 3, 0], [2, 1, 0, 3], [2, 1, 3, 0]]
    )
    assert asmnts[2][2:] == asmnts[3][2:]
