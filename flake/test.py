#!/usr/bin/env python3
#coding: utf-8
#                     The test module for the flake simulation
import unittest
from growth import Flake, Vector, VectorException, GridError

global F
F = Flake()


class TestLattice(unittest.TestCase):
  """ Flake with lattice, atoms and growth. """

  def test_grid(self):
    """ Test the raw_grid list and its access method grid().

    Creates new instance to ensure no altering of attributes. Ensure raw_grid and grid
    access method return False since no atom is present at initialization.
    """
    g = Flake(size=8)
    self.assertEqual(g.raw_grid[0][3][4][0], 0)
    self.assertEqual(g.grid(1, 2, 3), False)
    g.grid(3, 3, 3, 'set')      # Set atom
    self.assertEqual(g.grid(3, 3, 3), True)
    g.grid(3, 3, 3, 'remove')   # and remove it again
    self.assertEqual(g.grid(3, 3, 3), False)
    g.grid(3, 3, 3, 17)         # Set atom with specific energy parameter
    self.assertEqual(g.grid(3, 3, 3, 'get'), [True, 17])
    with self.assertRaises(GridError):
      g.grid(3, 3, 3, [2, 2])
    with self.assertRaises(GridError):
      g.grid(3, 3, 3, 'test_arg')


class TestVector(unittest.TestCase):

  def test_equality(self):
    vec = Vector(3, 4, 8)
    self.assertEqual(vec, Vector(3, 4, 8))

  def test_add(self):
    vec1 = Vector(3, 4, 7)
    vec2 = Vector(6, 3, 4)
    res = vec1 + vec2
    testvec = Vector(9, 7, 11)
    self.assertEqual(res, testvec)
    num = 17
    res = vec1 + num
    for (rcomp, vcomp) in zip(res, vec1):
      self.assertEqual(rcomp, vcomp + num)

  def test_sub(self):
    vec1 = Vector(2, 5, 0)
    vec2 = Vector(3, 2, -7)
    testvec = Vector(-1, 3, 7)
    res = vec1 - vec2
    self.assertEqual(res, testvec)
    num = 12
    res = vec1 - num
    for (rcomp, vcomp) in zip(res, vec1):
      self.assertEqual(rcomp, vcomp - num)

  def test_dist(self):
    vec1 = Vector(2, 5, 0)
    vec2 = Vector(3, 3, -7)
    self.assertEqual(7.3484692283495345, vec1.dist(vec2))
    self.assertEqual(vec1.dist(vec2), vec2.dist(vec1))
    self.assertEqual(0, vec1.dist(vec1))

  def test_dist_exception(self):
    with self.assertRaises(VectorException):
      Vector(1, 2, 3).dist("Failboat")

  def test_abs(self):
    vec1 = Vector(2, 5, 0)
    self.assertEqual(abs(vec1), vec1.dist(Vector(0, 0, 0)))


if __name__ == '__main__':
  unittest.main()
