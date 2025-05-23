"""
Test ReadScoring
"""

from whatshap.polyphase.solver import AlleleMatrix, scoreReadset
from whatshap.testhelpers import string_to_readset


def test_readscoring1():
    reads = """
    0001
     1001
      0101
       0100
        0110
         0001
          1001
    """
    am = AlleleMatrix(string_to_readset(reads), False)
    sim = scoreReadset(am, 2, 2)

    assert sim.get(0, 1) < 0.0
    assert sim.get(0, 2) > 0.0
    assert sim.get(0, 3) <= 0.0
    assert sim.get(0, 4) >= 0.0
    assert sim.get(0, 5) <= 0.0
    assert sim.get(0, 6) >= 0.0
    assert sim.get(1, 2) < 0.0
    assert sim.get(1, 3) > 0.0
    assert sim.get(1, 4) <= 0.0
    assert sim.get(1, 5) >= 0.0
    assert sim.get(1, 6) <= 0.0
    assert sim.get(2, 3) < 0.0
    assert sim.get(2, 4) > 0.0
    assert sim.get(2, 5) <= 0.0
    assert sim.get(2, 6) >= 0.0
    assert sim.get(3, 4) < 0.0
    assert sim.get(3, 5) > 0.0
    assert sim.get(3, 6) <= 0.0
    assert sim.get(4, 5) < 0.0
    assert sim.get(4, 6) > 0.0
    assert sim.get(5, 6) < 0.0


def test_readscoring2():
    reads = """
    0001
    0001
     1100
     1001
     11000
     11000
    """
    rs = string_to_readset(reads)
    am1 = AlleleMatrix(rs, False)
    am2 = AlleleMatrix(rs, True)
    sim1 = scoreReadset(am1, 2, 2)
    sim2 = scoreReadset(am2, 2, 2)

    assert len(am2) == 4

    assert sim1.get(0, 1) > 0.0
    assert sim1.get(0, 2) == sim1.get(1, 2) < 0.0
    assert sim1.get(0, 4) == sim1.get(1, 4) == sim1.get(0, 5) == sim1.get(1, 5) < 0.0
    assert sim2.get(0, 1) < 0.0
    assert sim2.get(0, 3) < 0.0
    assert sim2.get(0, 1) == 2 * sim1.get(0, 2)
    assert sim2.get(0, 3) == 4 * sim1.get(0, 5)
