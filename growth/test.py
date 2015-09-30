#!/usr/bin/env python3
#coding: utf-8
# The test module for our program

import unittest
from growth import *


class TestLattice(unittest.TestCase):
  """ Test the Flake. """

  def setUp(self):
    """ Instantiate the class. """
    self.Flake = Flake()

  def test_basis_vector(self):
    """ Base Vector of the lattice. """
    self.assertTrue(len(self.Flake.fcc_base) == 3)
    a = (1, 0, 0)
    b = (1/sqrt(2), 1/sqrt(2), 0)
    c = (1/2, 1/(2*sqrt(3)), sqrt(2/3))
    test_vecs = np.array((a, b, c)).T
    f_base_vec = self.Flake.fcc_base
    for (idx, vec) in enumerate(test_vecs):
      self.assertTrue(np.array_equal(f_base_vec[idx], vec))

  def test_constructor(self):
    """ Check if the constructor creates a vector correctly and can reach all
    lattice points.
    """
    # Testing Varibales:
    steps = np.array((1, 10, 0))
    test_base = np.array(((1, 2, 3), (1, 0, 1), (2, 2, 1)))
    result_1 = self.Flake.constructor(steps, test_base)
    exp = np.array((21, 1, 22)).T
    # print("exp, res:\n", exp, result_1)
    self.assertTrue(np.array_equal(result_1, np.array(exp)))

    int_coords = (17, 8, -3)
    x = (17 + 8*(1/sqrt(2)) + 0.5*-3)
    y = (8/sqrt(2) - 3/(2*sqrt(3)))
    z = (-3*sqrt(2/3))
    expected = np.array((x, y, z))
    test_result2 = self.Flake.constructor(int_coords, basis='full')
    self.assertTrue(np.allclose(np.array(expected).T, test_result2))

  def test_plot(self):
    """ Very basic testing of the plot method. First plot should return points
    on a diagonal.
    """
    # return  # THIS ALWAYS PASSES! Just to not plot all the time
    testpoints = np.array(((1, 1, 1), (2, 2, 2), (3, 3, 3), (4, 4, 4)))
    self.assertIsNone(self.Flake.plot(testpoints))
    self.assertIsNone(self.Flake.plot())

  def test_permutator(self):
    input = (1, 0)
    output = self.Flake.permutator(input)
    expected = ((1, 1), (0, 1), (1, 0), (0, 0))
    self.assertEqual(output, expected)

    input = range(-1, 1)
    output = self.Flake.permutator(input, inplane=False)
    expected = ((-1, -1, -1), (0, -1, -1), (-1, 0, -1), (-1, -1, 0),
                (-1, 0, 0), (0, 0, -1), (0, -1, 0), (0, 0, 0))
    self.assertEqual(output, expected)

  # def test_stack(self):
    # output = self.Flake.stack(height=5, size=3)
    # self.Flake.plot(output)

  # def test_next_neighbour(self):
    # return
    # inplane = [(1, 0, 0),  (-1, 0, 0), (0, 1, 0),  (0, -1, 0),
               # (1, -1, 0), (-1, 1, 0)]
    # f2f = inplane.extend([(1, 0, -1), (-1, 0, 1), (0, 1, -1), (0, -1, 1),
                          # (0, 0, 1),  (0, 0, -1)])

    # f2h = inplane.extend([(0, 1, -1), (0, -1, 1), (0, 0, 1),  (0, 0, -1),
                          # (1, -1, 1), (1, 0, -1)])
    # self.assertEqual(self.Flake.next_neighbour('f2f'), f2f)
    # self.assertEqual(self.Flake.next_neighbour('h2h'), h2h)
    # self.assertEqual(self.Flake.next_neighbour('f2h'), f2h)
    # self.assertEqual(self.Flake.next_neighbour('h2f'), h2f)

if __name__ == '__main__':
  unittest.main()
