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
    res1 = func(*index, twin=0)
    exp1 = self.F.Vector(3, 3.4641016151377544, 4.898979485566356)
    res2 = func(*index, twin=1).y
    exp2 = sqrt(3)*(2+2 * 1/3)
    self.assertEqual(res1, exp1)
    self.assertEqual(res2, exp2)

  def test_layer_generator(self):
    """ Sets the configuration of twin planes. """
    g = Flake(14)
    twins = (2, 3, 5)
    layers = g.layer_generator(twins)
    self.assertTrue(all(layers[i] == 0 for i in range(0, 3)))
    self.assertTrue(all(layers[i] == 1 for i in range(3, 4)))
    self.assertTrue(all(layers[i] == 2 for i in range(4, 6)))
    self.assertTrue(all(layers[i] == 3 for i in range(7, 10)))

  def test_nn_gen(self):
    """ Next neighbours vector generator. """
    pass

if __name__ == '__main__':
  unittest.main()
