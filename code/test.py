#!/usr/bin/env python3
#coding: utf-8
# The test module for our program
# pylint: disable=W0401

import unittest
from growth import *


class TestLattice(unittest.TestCase):
  """ Test the Flake. """

  def setUp(self):
    """ Instantiate the class. """
    self.F = Flake()

  def test_basis_vector(self):
    """ Base Vector of the lattice. """
    a = (1, 0, 0)
    b = (1/sqrt(2), 1/sqrt(2), 0)
    c = (1/2, 1/(2*sqrt(3)), sqrt(2/3))
    expected = self.F.LatticeBase(a, b, c)
    output = self.F.fcc_base
    self.assertEqual(expected, output)
    self.assertTrue(len(self.F.fcc_base) == 3)

  def test_grid(self):
    """ Test the boolean grid_list of atoms. """
    g = Flake()   # Create new Flake instance of to ensure no altering of inner
                  # values.
    self.assertEqual(g.grid_list[0][3][4], 0)
    self.assertEqual(g.grid(1, 2, 3, None), False)
    # self.assertRaises(TypeError, g.grid(3, 3, 3, 'that_string'))
    grid_point = g.grid(4, 1, 2, val=1)
    print("GP", grid_point, type(grid_point))
    self.assertTrue(type(grid_point) is bool)
    index = (2,1,4)
    g.grid(*index, val=1)
    self.assertEqual(g.grid(*index), True)

if __name__ == '__main__':
  unittest.main()
