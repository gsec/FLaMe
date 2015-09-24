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
    self.assertTrue(len(self.Flake.basis_vector) == 3)
    a = (1, 0, 0)
    b = (1/sqrt(2), 1/sqrt(2), 0)
    c = (1/2, 1/(2*sqrt(3)), sqrt(2/3))

    test_vecs = self.Flake.Vector(a, b, c)
    f_base_vec = self.Flake.basis_vector

    for (idx, vec) in enumerate(test_vecs):
      self.assertTrue(f_base_vec[idx] == vec)

if __name__ == '__main__':
  unittest.main()
