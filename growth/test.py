#!/usr/bin/env python3
#coding: utf-8
# The test module for our program
# pylint: disable=W0401

import unittest
from growth2 import *


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
    expected = self.F.LatticeBase(a=a, b=b, c=c)
    output = self.F.fcc_base
    self.assertEqual(expected, output)
    self.assertTrue(len(self.F.fcc_base) == 3)

  def test_grid(self):
    """ Test the boolean grid of atoms. """
    self.assertEqual(self.F.grid_list[0][3][4], 0)
    index = (2,1,4)
    self.F.grid(*index, val=1)
    self.assertEqual(self.F.grid(*index), True)

if __name__ == '__main__':
  unittest.main()
