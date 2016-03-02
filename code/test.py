#!/usr/bin/env python
#coding: utf-8
#                     The test module for the flake simulation
from __future__ import print_function, division, generators
import unittest
import growth as g
from helper import Vector


class TestFlake(unittest.TestCase):
  """ Flake sanity checks.
  """
  def test_point_structure(self):
    f = g.Flake(seed='point')
    fn = set([(-1, 0, 1), (0, 0, 1), (0, -1, 1),
              (0, -1, 0), (-1, 0, 0), (0, 1, 0),
              (-1, 1, 0), (1, 0, 0), (-1, -1, 0),
              (0, -1, -1), (0, 0, -1), (-1, -1, -1)])
    self.assertEqual(f.atoms, set([(0, 0, 0)]))
    self.assertEqual(f.surface[1], fn)
    self.assertEquals(len(f.integrated_surface()), len(f.maxNB))


class TestVector(unittest.TestCase):
  def test_equality(self):
    vec = Vector((3, 4, 8))
    self.assertEqual(vec, Vector((3, 4, 8)))

  def test_add(self):
    vec1 = Vector((3, 4, 7))
    vec2 = Vector((6, 3, 4))
    res = vec1 + vec2
    testvec = Vector((9, 7, 11))
    self.assertEqual(res, testvec)

  def test_sub(self):
    vec1 = Vector((2, 5, 0))
    vec2 = Vector((3, 2, -7))
    testvec = Vector((-1, 3, 7))
    res = vec1 - vec2
    self.assertEqual(res, testvec)

  def test_dist(self):
    vec1 = Vector((2, 5, 0))
    vec2 = Vector((3, 3, -7))
    self.assertAlmostEqual(7.3484692283495345, vec1.dist(vec2))
    self.assertEqual(vec1.dist(vec2), vec2.dist(vec1))
    self.assertEqual(0, vec1.dist(vec1))

  def test_abs(self):
    vec1 = Vector((2, 5, 0))
    self.assertEqual(abs(vec1), vec1.dist(Vector((0, 0, 0))))


if __name__ == '__main__':
  unittest.main()
