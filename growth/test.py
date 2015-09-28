#!/usr/bin/env python3
#coding: utf-8
# The test module for our program

import unittest
from growth import *


class TestLattice(unittest.TestCase):
  """ the first test."""

  def setUp(self):
    """ Instantiate the class."""
    self.Flake = Flake()

  def test_basis_vector(self):
    """ Base Vector of the lattice.
    """
    self.assertTrue(len(self.Flake.basis_vector) == 3)
    a = (1, 0, 0)
    b = (1/sqrt(2), 1/sqrt(2), 0)
    c = (1/2, 1/(2*sqrt(3)), sqrt(2/3))

    test_vecs = matrix((a, b, c))
    f_base_vec = self.Flake.basis_vector

    for (idx, vec) in enumerate(test_vecs):
      self.assertTrue(np.array_equal(f_base_vec[idx], vec))

  def test_constructor(self):
    """ Check if the constructor creates a vector correctly and can reach all
    lattice points.
    """
    ordinal = matrix((1, 10, 0))
    test_set = matrix(((1, 2, 3), (1, 0, 1), (2, 2, 1)))
    test_result1 = self.Flake.coordinates(ordinal, test_set)
    self.assertTrue(np.array_equal(test_result1, matrix((11, 2, 13))))
    int_coords = (17, 8, -3)
    x = (17 + 8*(1/sqrt(2)) + 0.5*-3)
    y = (8/sqrt(2) - 3/(2*sqrt(3)))
    z = (-3*sqrt(2/3))
    expected = matrix((x, y, z))
    test_result2 = self.Flake.coordinates(int_coords)
    self.assertTrue(np.allclose(expected, test_result2))

  def test_plot(self):
    """ Very basic testing of the plot method. First plot should return points
    on a diagonal.
    """
    testpoints = matrix(((1, 1, 1), (2, 2, 2), (3, 3, 3), (4, 4, 4)))
    # THIS ALWAYS PASSES! Just to not plot all the time
    # return True
    self.assertIsNone(self.Flake.plot(testpoints))
    self.assertIsNone(self.Flake.plot())

  def test_plane(self):
    size = 4
    z = 0
    layer = self.Flake.plane(z=z, size=size)
    print(layer)

if __name__ == '__main__':
  unittest.main()
