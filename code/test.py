#!/usr/bin/env python
#coding: utf-8
#                     The test module for the flake simulation
from __future__ import print_function, division, generators
import unittest
import growth
from helper import Vector


class TestFlake(unittest.TestCase):
  """ Flake sanity checks.
  """
  def test_point_structure(self):
    f = growth.Flake(seed='point')
    surface_one = set([(-1, 0, 1), (0, 0, 1), (0, -1, 1),
              (0, -1, 0), (-1, 0, 0), (0, 1, 0),
              (-1, 1, 0), (1, 0, 0), (-1, -1, 0),
              (0, -1, -1), (0, 0, -1), (-1, -1, -1)])
    self.assertEqual(f.atoms, set([(0, 0, 0)]))
    self.assertEqual(f.surface[1], surface_one)
    self.assertEquals(len(f.integrated_surface()), len(f.maxNB))


  def test_neighbours(self):
    """ Test if our chosen atom has the proper neighbours in index space,
    ensure this translates into correct atoms considering the real space
    distance and make sure all neighbours return our seed as non-void
    neighborhood.
    """
    f = growth.Flake(seed='point')
    atom = (2, 5, 8)
    all_neighbours = [(2, 6, 9), (3, 4, 7), (3, 5, 9), (3, 4, 8), (1, 6, 7),
                      (3, 6, 7), (2, 4, 7), (2, 5, 7), (1, 5, 9), (1, 5, 8),
                      (3, 5, 8), (1, 5, 7), (2, 5, 9), (2, 4, 9), (2, 4, 8),
                      (3, 6, 9), (1, 6, 8), (3, 6, 8), (1, 6, 9), (1, 4, 9),
                      (3, 5, 7), (1, 4, 8), (3, 4, 9), (2, 6, 8), (2, 6, 7),
                      (1, 4, 7)]
    self.assertEqual(all_neighbours, f.abs_neighbours(atom))

    point_seed = (0, 0, 0)
    next_neighbours = [(0, -1, -1), (0, 0, -1),  (-1, 0, 1),  (-1, 0, 0),
                       (1, 0, 0),   (0, 0, 1),   (0, -1, 1),  (0, -1, 0),
                       (-1, 1, 0),  (-1, -1, 0), (0, 1, 0),   (-1, -1, -1)]
    self.assertEqual(next_neighbours, f.real_neighbours(point_seed, void=True))
    for each in next_neighbours:
      self.assertEqual([point_seed], f.real_neighbours(each, void=False))


  def test_twin_plane_creation(self):
    twinplanes = (4, 5, 6)
    self.assertEqual(twinplanes, (000))


class TestVector(unittest.TestCase):
  """ Sanity checks of the Vector class for basic operations.
  """
  def test_equality(self):
    vec = Vector(3, 4, 8)
    self.assertEqual(vec, Vector(3, 4, 8))

  def test_add(self):
    vec1 = Vector(3, 4, 7)
    vec2 = Vector(6, 3, 4)
    res = vec1 + vec2
    testvec = Vector(9, 7, 11)
    self.assertEqual(res, testvec)

  def test_sub(self):
    vec1 = Vector(2, 5, 0)
    vec2 = Vector(3, 2, -7)
    testvec = Vector(-1, 3, 7)
    res = vec1 - vec2
    self.assertEqual(res, testvec)

  def test_dist(self):
    vec1 = Vector(2, 5, 0)
    vec2 = Vector(3, 3, -7)
    self.assertEqual(7.3484692283495345, vec1.dist(vec2))
    self.assertEqual(vec1.dist(vec2), vec2.dist(vec1))
    self.assertEqual(0, vec1.dist(vec1))

  def test_abs(self):
    vec1 = Vector(2, 5, 0)
    self.assertEqual(abs(vec1), vec1.dist(Vector(0, 0, 0)))


if __name__ == '__main__':
  unittest.main()
