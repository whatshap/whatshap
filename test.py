#!/usr/bin/env python3

from whatshap._core import PyRead, PyReadSet

s = PyReadSet()

r = PyRead('Read A', 56)
r.addVariant(100,'A', 1, 37)
r.addVariant(101,'C', 0, 18)
s.add(r)

r = PyRead('Read B', 0)
r.addVariant(100,'C', 0, 23)
s.add(r)

print(s)
