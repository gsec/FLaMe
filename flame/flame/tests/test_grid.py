#!/usr/bin/env python
#coding: utf-8

""" Tests for the `Grid` module. Check for basic math and sanity of our low-level module.
"""

from __future__ import print_function, division, generators
import unittest
from flame.grid import Grid, Vector

ATOM_DIA = 2


class TestPerfectGrid(unittest.TestCase):
    """ Grid base tests on what should be expected for simple configurations.
    """
    def setUp(self):
        self.pGrid = Grid(())           # create a perfect Grid
        self.zero_vector = Vector(0, 0, 0)

    def test_twin_gen(self):
        self.assertEqual(self.pGrid.twin_gen(), (None, None, None))

        for i in range(-20, 20):
            self.assertEqual(self.pGrid.shift(i), i % 3)

    def test_coords(self):
        center_atom = self.pGrid.coord((0, 0, 0))
        self.assertEqual(center_atom, self.zero_vector)
        atom_set = {(-1, 0, 0),
                    (0, -1, 0),
                    (0, 0, -1),
                    (1, 0, 0),
                    (0, 1, 0),
                    (0, 0, 1)}
        atom_coords = (self.pGrid.coord(at) for at in atom_set)
        for atom in atom_coords:
            self.assertAlmostEqual(atom.dist(center_atom), ATOM_DIA)


class TestTwinnedGrid(unittest.TestCase):
    """ Grid base tests on what should be expected for simple configurations.
    """
    def setUp(self):
        self.tGrid = Grid((-1, 1))      # create the flawed Grid with twin planes
        self.zero_vector = Vector(0, 0, 0)

    def test_twin_gen(self):
        sign = 1
        shift = 0
        self.assertEqual(self.tGrid.twin_gen(), ([2, 1, 0], shift, sign))

        for i in range(-20, -1):
            self.assertEqual(self.tGrid.shift(i), i % 3)


if __name__ == '__main__':
    unittest.main()
