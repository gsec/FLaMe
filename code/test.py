#!/usr/bin/env python3
#coding: utf-8
# The test module for our program

import unittest
from growth import Flake, sqrt


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

    test_vecs = self.Flake.Vector((a, b, c))
    f_base_vec = self.Flake.basis_vector

    for (idx, vec) in enumerate(test_vecs):
      # print("YOLO:", f_base_vec[idx]-vec)
      # print("IDX:", f_base_vec[idx], "\t\t", "vec:", vec)
      self.assertEqual(f_base_vec[idx].all(), vec.all())

  def test_constructor(self):
    """ Check if the constructor creates a vector correctly and can reach all
    lattice points.
    """
    int_coords = (17, 8, -3)
    x = (17 + 8*(1/sqrt(2)) + 0.5*-3)
    y = (8/sqrt(2) - 3/(2*sqrt(2)))
    z = (-3*sqrt(2/3))
    expected = self.Flake.Vector((x, y, z))
    test_result = self.Flake.coordinates(int_coords)
    # print("Exp: ", expected, "\nResult:", test_result)
    for (e, t) in zip(test_result, expected):
      print(e, t)
      self.assertEqual(e, t)
    # (self.assertEqual(e, t) for (e, t) in zip(test_result, expected))

if __name__ == '__main__':
  unittest.main()
