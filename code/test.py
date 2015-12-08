#!/usr/bin/env python
#coding: utf-8
#                     The test module for the flake simulation
from __future__ import print_function, division, generators
import unittest
from helper import Vector
from np import Grid


class TestVector(unittest.TestCase):
  def test_equality(self):
    vec = Vector((3, 4, 8))
    print(vec)
    self.assertEqual(vec, Vector((3, 4, 8)))

  def test_add(self):
    vec1 = Vector((3, 4, 7))
    vec2 = Vector((6, 3, 4))
    res = vec1 + vec2
    testvec = Vector((9, 7, 11))
    self.assertEqual(res, testvec)
    num = 17
    res = vec1 + num
    for (rcomp, vcomp) in zip(res, vec1):
      self.assertEqual(rcomp, vcomp + num)

  def test_sub(self):
    vec1 = Vector((2, 5, 0))
    vec2 = Vector((3, 2, -7))
    testvec = Vector((-1, 3, 7))
    res = vec1 - vec2
    self.assertEqual(res, testvec)
    num = 12
    res = vec1 - num
    for (rcomp, vcomp) in zip(res, vec1):
      self.assertEqual(rcomp, vcomp - num)

  def test_dist(self):
    vec1 = Vector((2, 5, 0))
    vec2 = Vector((3, 3, -7))
    self.assertEqual(7.3484692283495345, vec1.dist(vec2))
    self.assertEqual(vec1.dist(vec2), vec2.dist(vec1))
    self.assertEqual(0, vec1.dist(vec1))

  def test_dist_exception(self):
    with self.assertRaises(TypeError):
      Vector((1, 2, 3)).dist("Failboat")

  def test_abs(self):
    vec1 = Vector((2, 5, 0))
    self.assertEqual(abs(vec1), vec1.dist(Vector((0, 0, 0))))


class SimpleFlakeTest(unittest.TestCase):
  def setUp(self):
    self.G = Grid(size=20, twins=(3,))


class TestGrid(SimpleFlakeTest):
  def test_init(self):
    conditions = ((len(self.G.core[0]), 20), (self.G.core.ndim, 4),
    (self.G.size, 20), (len(self.G.core[2, 3, 4]), 4))
    for res, exp in conditions:
      self.assertEqual(res, exp)

  def test_seed(self):
    seeds = self.G.permutator((9, 10, 11))
    for i, j, k in seeds:
      self.assertEqual(self.G.core[i, j, k, 0], 0)


if __name__ == '__main__':
  unittest.main()
