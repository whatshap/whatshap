"""
Test ReadScoring
"""

from whatshap.core import Read, ReadSet
from whatshap.polyphase.solver import AlleleMatrix, scoreReadset


def test_readscoring_toy():
    readset = ReadSet()
    read1 = Read("name1", 15)
    read1.add_variant(0, 0, 1)
    read1.add_variant(1, 0, 1)
    read1.add_variant(2, 0, 1)
    read1.add_variant(3, 1, 1)
    readset.add(read1)
    read2 = Read("name2", 15)
    read2.add_variant(1, 1, 1)
    read2.add_variant(2, 0, 1)
    read2.add_variant(3, 0, 1)
    read2.add_variant(4, 1, 1)
    readset.add(read2)
    read3 = Read("name3", 15)
    read3.add_variant(2, 0, 1)
    read3.add_variant(3, 1, 1)
    read3.add_variant(4, 0, 1)
    read3.add_variant(5, 1, 1)
    readset.add(read3)
    read4 = Read("name4", 15)
    read4.add_variant(3, 0, 1)
    read4.add_variant(4, 1, 1)
    read4.add_variant(5, 0, 1)
    read4.add_variant(6, 0, 1)
    readset.add(read4)
    read5 = Read("name5", 15)
    read5.add_variant(4, 0, 1)
    read5.add_variant(5, 1, 1)
    read5.add_variant(6, 1, 1)
    read5.add_variant(7, 0, 1)
    readset.add(read5)
    read6 = Read("name6", 15)
    read6.add_variant(5, 0, 1)
    read6.add_variant(6, 0, 1)
    read6.add_variant(7, 0, 1)
    read6.add_variant(8, 1, 1)
    readset.add(read6)
    read7 = Read("name7", 15)
    read7.add_variant(6, 1, 1)
    read7.add_variant(7, 0, 1)
    read7.add_variant(8, 0, 1)
    read7.add_variant(9, 1, 1)
    readset.add(read7)
    am = AlleleMatrix(readset)
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
