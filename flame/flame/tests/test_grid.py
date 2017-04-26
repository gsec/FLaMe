""" Tests for the `Grid` module. Check for basic math and sanity of our low-level module.
"""
from __future__ import print_function, division, generators
import unittest
from flame.grid import Grid, Vector

ATOM_DIA = 2


class TestVectorMath(unittest.TestCase):
    def setUp(self):
        self.vec0 = Vector(2, 3, 4)
        self.vec1 = Vector(3, 4, 8)

    def test_add(self):
        vec1 = Vector(3, 4, 7)
        vec2 = Vector(6, 3, 4)
        res = vec1 + vec2
        testvec = Vector(9, 7, 11)
        self.assertEqual(res, testvec)

    def test_basic(self):
        self.assertEqual("Vector:(2.0, 3.0, 4.0)", str(self.vec0))
        self.assertEqual(self.vec0, Vector(2, 3, 4))
        for comp, num in zip(self.vec0, (2, 3, 4)):
            self.assertEqual(comp, num)

    def test_equality(self):
        self.assertEqual(self.vec1, Vector(3, 4, 8))
        self.assertNotEqual(self.vec1, Vector(1, 2, 3))

    def test_scaling(self):
        self.assertEqual(self.vec0/2, Vector(1, 1.5, 2))
        self.assertEqual(self.vec0*2, Vector(4, 6, 8))
        self.assertEqual(self.vec1*2, 2*self.vec1)
        with self.assertRaises(AssertionError):
            self.vec0*self.vec0
        with self.assertRaises(AssertionError):
            self.vec0/Vector(1, 2, 3)

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
