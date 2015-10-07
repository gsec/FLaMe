#!/usr/bin/env python3
#coding: utf-8
# The test module for our program
# pylint: disable=W0401

import unittest
from growth import *


class TestLattice(unittest.TestCase):
  """ Flake with lattice, atoms and growth. """

  def setUp(self):
    """ Instantiate the class. """
    self.F = Flake()

  def test_basis_vector(self):
    """ FCC base vectors of the lattice (a,b,c) """
    a = (1, 0, 0)
    b = (1/sqrt(2), 1/sqrt(2), 0)
    c = (1/2, 1/(2*sqrt(3)), sqrt(2/3))
    expected = self.F.LatticeBase(a, b, c)
    output = self.F.fcc_base
    self.assertEqual(expected, output)
    self.assertTrue(len(self.F.fcc_base) == 3)

  def test_grid(self):
    """ Boolean n*n*n list of atoms.
    Creates new instance to ensure no altering of attributes."""
    g = Flake()
    self.assertEqual(g.grid_list[0][3][4], 0)
    self.assertEqual(g.grid(1, 2, 3, None), False)
    self.assertEqual(g.grid(0, 2, 0), False)

    grid_point = g.grid(4, 1, 2, val=1)
    # print("GP", grid_point, type(grid_point))
    self.assertTrue(isinstance(grid_point, bool))
    index = (2, 1, 4)
    g.grid(*index, val=1)
    self.assertEqual(g.grid(*index), True)
    with self.assertRaises(TypeError):
      g.grid(3, 3, 3, 'that_string')

  def test_coord(self):
    """ Retrieve cartesian coordinates from lattice indices.
    """
    func = self.F.coord
    index = (1, 2, 3)
    self.assertEqual(func(*index), self.F.Vector(3, 4.041451884327381,
                                              4.898979485566356))
    self.assertEqual(func(*index, struct='hcp'),
                     self.F.Vector(3, 3.4641016151377544, 4.898979485566356))
    self.assertEqual(func(*index), func(*index, struct='fcc'))
    with self.assertRaises(KeyError):
      func(2, 3, 4, 'nonexistent_structure')

  def test_displacement_vector(self):
    """ Displacement vector between `fcc` and `hcp` plane.
    """
    delta_expected = -0.5773502691896262
    delta_output = self.F.coord(2, 2, 3, struct='hcp').y - self.F.coord(
      2, 2, 3, 'fcc').y
    expected_lattice_delta = self.F.Vector(0, delta_expected, 0)
    self.assertEqual(delta_expected, delta_output)
    self.assertEqual(expected_lattice_delta, self.F.lattice_delta)

if __name__ == '__main__':
  unittest.main()
