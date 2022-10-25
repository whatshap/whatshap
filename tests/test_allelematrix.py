from whatshap.polyphase.solver import AlleleMatrix
from whatshap.testhelpers import string_to_readset


def test_get_allele1():
    reads = """
    001001
    110101
    """
    am = AlleleMatrix(string_to_readset(reads))
    assert len(am) == 2
    assert am.getNumPositions() == 6
    assert am.getAllele(0, 0) == 0
    assert am.getAllele(0, 2) == 1
    assert am.getAllele(0, 6) == -1
    assert am.getAllele(1, 0) == 1
    assert am.getAllele(1, 6) == -1


def test_get_allele2():
    reads = """
    110101
     01  01
      001001
    """
    am = AlleleMatrix(string_to_readset(reads))
    assert len(am) == 3
    assert am.getNumPositions() == 8
    assert am.getAllele(2, 0) == -1
    assert am.getAllele(2, 2) == 0
    assert am.getAllele(2, 9) == -1
    assert am.getAllele(0, 0) == 1
    assert am.getAllele(0, 2) == 0
    assert am.getAllele(0, 3) == 1
    assert am.getAllele(0, 7) == -1
    assert am.getAllele(1, 0) == -1
    assert am.getAllele(1, 2) == 1
    assert am.getAllele(1, 3) == -1
    assert am.getAllele(1, 6) == 1
    assert am.getAllele(1, 7) == -1


def test_get_read1():
    reads = """
    110101
     01  01
      001001
    """
    am = AlleleMatrix(string_to_readset(reads))
    assert am.getRead(0) == [(0, 1), (1, 1), (2, 0), (3, 1), (4, 0), (5, 1)]
    assert am.getRead(1) == [(1, 0), (2, 1), (5, 0), (6, 1)]
    assert am.getRead(2) == [(2, 0), (3, 0), (4, 1), (5, 0), (6, 0), (7, 1)]
    assert am.getFirstPos(0) == 0 < 5 == am.getLastPos(0)
    assert am.getFirstPos(1) == 1 < 6 == am.getLastPos(1)
    assert am.getFirstPos(2) == 2 < 7 == am.getLastPos(2)


def test_get_positions1():
    reads = """
    1 101 01
      01   01
       00 1001
    """
    am = AlleleMatrix(string_to_readset(reads))
    gen_positions = [10 * (x + 1) for x in [0, 2, 3, 4, 6, 7, 8, 9]]
    assert am.getPositions() == gen_positions
    for pos, genpos in enumerate(gen_positions):
        assert am.globalToLocal(genpos) == pos
        assert am.localToGlobal(pos) == genpos


def test_get_alleledepths1():
    reads = """
    1 101 01
      01   01
       00 1001
        1 00 1
    """
    am = AlleleMatrix(string_to_readset(reads))
    ad = [[0, 1], [1, 1], [2, 1], [1, 2], [2, 1], [3, 1], [1, 1], [0, 2]]
    for i in range(am.getNumPositions()):
        assert am.getAlleleDepths(i) == ad[i]


def test_get_alleledepths2():
    reads = """
    1 101 01
      01   01
       00 1001
        1 00 2
    """
    am = AlleleMatrix(string_to_readset(reads))
    ad = [[0, 1, 0], [1, 1, 0], [2, 1, 0], [1, 2, 0], [2, 1, 0], [3, 1, 0], [1, 1, 0], [0, 1, 1]]
    for i in range(am.getNumPositions()):
        assert am.getAlleleDepths(i) == ad[i]


def test_get_alleledepths3():
    reads = """
    1 101 01
      01   01
       00 1001
        1 00 2
    """
    am = AlleleMatrix(string_to_readset(reads))
    ad = [[0, 1, 0], [1, 1, 0], [2, 1, 0], [1, 2, 0], [2, 1, 0], [3, 1, 0], [1, 1, 0], [0, 1, 1]]
    for i in range(am.getNumPositions()):
        assert am.getAlleleDepths(i) == ad[i]


def test_sub_interval1():
    reads = """
    1001 01001
      1010010  01
       100 10 0010
        010  100
          001 100 01
    """
    am = AlleleMatrix(string_to_readset(reads))
    s1 = am.extractInterval(0, 16)
    assert len(am) == len(s1)
    assert am.getNumPositions() == s1.getNumPositions()
    for i in range(len(am)):
        for j in range(am.getNumPositions()):
            assert am.getAllele(i, j) == s1.getAllele(i, j)
            assert am.getAlleleDepths(j) == s1.getAlleleDepths(j)
            assert am.localToGlobal(j) == s1.localToGlobal(j)

    s2 = am.extractInterval(2, 13)
    assert len(s2) == 5
    assert s2.getNumPositions() == 11
    for i in range(len(s2)):
        for j in range(s2.getNumPositions()):
            assert am.getAllele(i, j + 2) == s2.getAllele(i, j)
            assert am.getAlleleDepths(j + 2) == s2.getAlleleDepths(j)
            assert am.localToGlobal(j + 2) == s2.localToGlobal(j)
    assert s2.getRead(0) == [(0, 0), (1, 1), (3, 0), (4, 1), (5, 0), (6, 0), (7, 1)]


def test_sub_interval2():
    reads = """
    1001 01001
      1010010  01
       100 10 0010
        010  100
          001 100 01
    """
    am = AlleleMatrix(string_to_readset(reads))
    s1 = am.extractInterval(0, 4, True)
    s2 = am.extractInterval(0, 4, False)
    s3 = am.extractInterval(10, 16, True)
    s4 = am.extractInterval(10, 16, False)
    assert len(am) == len(s2) == len(s4)
    assert len(s1) == 3
    assert len(s3) == 4
    assert s1.getNumPositions() == s2.getNumPositions() == 4
    assert s3.getNumPositions() == s4.getNumPositions() == 6
    for i in range(len(s3)):
        for j in range(s3.getNumPositions()):
            assert am.getAllele(i + 1, j + 10) == s3.getAllele(i, j)
            assert am.getAlleleDepths(j + 10) == s3.getAlleleDepths(j)
            assert am.localToGlobal(j + 10) == s3.localToGlobal(j)
    assert s4.getRead(0) == []
    assert s2.getRead(3) == s2.getRead(4) == []


def test_sub_matrix1():
    reads = """
    1001 01001
      1010010  01
       100 10 0010
        010  100
          001 100 01
    """
    am = AlleleMatrix(string_to_readset(reads))
    pos = [0, 1, 9, 10, 13, 14]
    s1 = am.extractSubMatrix(pos, [0, 1, 2, 3, 4], True)
    s2 = am.extractSubMatrix(pos, [0, 1, 2, 3, 4], False)
    assert len(am) == len(s2)
    assert len(s1) == 4
    assert s1.getNumPositions() == s2.getNumPositions() == 6
    for i in range(len(s2)):
        for j in range(s2.getNumPositions()):
            assert am.getAllele(i, pos[j]) == s2.getAllele(i, j)
            assert am.getAlleleDepths(pos[j]) == s2.getAlleleDepths(j)
            assert am.localToGlobal(pos[j]) == s2.localToGlobal(j)
    assert s1.getRead(0) == s2.getRead(0)
    assert s1.getRead(1) == s2.getRead(2)
    assert s1.getRead(2) == s2.getRead(3)
    assert s1.getRead(3) == s2.getRead(4)


def test_sub_matrix2():
    reads = """
    1001 01001
      1010010  01
       100 10 0010
        010  100
          001 100 01
    """
    am = AlleleMatrix(string_to_readset(reads))
    pos1 = list(range(16))
    pos2 = [0, 1, 9, 10, 13, 14]
    reads = [1, 2, 3]
    s1 = am.extractSubMatrix(pos1, reads, True)
    s2 = am.extractSubMatrix(pos2, reads, True)
    assert len(s1) == 3
    assert len(s2) == 2
    assert s1.getNumPositions() == 16
    assert s2.getNumPositions() == 6
    for i in range(len(s1)):
        for j in range(s1.getNumPositions()):
            assert am.getAllele(i + 1, j) == s1.getAllele(i, j)
    for i in range(len(s2)):
        for j in range(s2.getNumPositions()):
            assert am.getAllele(i + 2, pos2[j]) == s2.getAllele(i, j)
