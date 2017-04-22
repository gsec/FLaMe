#coding: utf-8

import unittest
from flame import growth
from flame.grid import Vector


class TestFlakeBasics(unittest.TestCase):
  """ Flake sanity checks.
  """
  def test_seed_init(self):
    f = growth.Flake(seed='point')
    surface_one = set([(-1, 0, 1), (0, 0, 1), (0, -1, 1), (0, -1, 0),
                        (-1, 0, 0), (0, 1, 0), (-1, 1, 0), (1, 0, 0),
                        (-1, -1, 0), (0, -1, -1), (0, 0, -1), (-1, -1, -1)])
    self.assertEqual(f.atoms, set([(0, 0, 0)]))
    self.assertEqual(f.surface[1], surface_one)
    self.assertEqual(len(f.integrated_surface()), len(f.maxNB))

    self.assertEqual(f.grid.twins, ())
    self.assertEqual(f.grid.twin_layers, None)


  def test_neighbours(self):
    """ Test if our chosen atom has the proper neighbours in index space,
    ensure this translates into correct atoms considering the real space
    distance and make sure all neighbours return our seed as non-void
    neighborhood.
    """
    f = growth.Flake(seed='point')
    atom = (2, 5, 8)
    all_neighbours = [(2, 6, 9), (3, 4, 7), (3, 5, 9), (3, 4, 8), (1, 6, 7),
                      (3, 6, 7), (2, 4, 7), (2, 5, 7), (1, 5, 9), (1, 5, 8),
                      (3, 5, 8), (1, 5, 7), (2, 5, 9), (2, 4, 9), (2, 4, 8),
                      (3, 6, 9), (1, 6, 8), (3, 6, 8), (1, 6, 9), (1, 4, 9),
                      (3, 5, 7), (1, 4, 8), (3, 4, 9), (2, 6, 8), (2, 6, 7),
                      (1, 4, 7)]
    self.assertEqual(set(all_neighbours), set(f.abs_neighbours(atom)))

    point_seed = (0, 0, 0)
    next_neighbours = [(0, -1, -1), (0, 0, -1),  (-1, 0, 1),  (-1, 0, 0),
                       (1, 0, 0),   (0, 0, 1),   (0, -1, 1),  (0, -1, 0),
                       (-1, 1, 0),  (-1, -1, 0), (0, 1, 0),   (-1, -1, -1)]
    self.assertEqual(set(next_neighbours), set(f.real_neighbours(point_seed, void=True)))
    for each in next_neighbours:
      self.assertEqual([point_seed], f.real_neighbours(each, void=False))


  def test_twin_plane_creation(self):
    twinplanes = (4, 5, 6)
    f = growth.Flake(*twinplanes)
    self.assertEqual(f.grid.twins, twinplanes)

    layers = [1, 0, 1]
    self.assertEqual(f.grid.twin_layers, layers)


class TestVector(unittest.TestCase):
  """ Sanity checks of the Vector class for basic operations.
  """
  def test_equality(self):
    vec = Vector(3, 4, 8)
    self.assertEqual(vec, Vector(3, 4, 8))
    self.assertNotEqual(vec, Vector(1, 2, 3))

  def test_add(self):
    vec1 = Vector(3, 4, 7)
    vec2 = Vector(6, 3, 4)
    res = vec1 + vec2
    testvec = Vector(9, 7, 11)
    self.assertEqual(res, testvec)

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


class TestFlakeGrowth(unittest.TestCase):
    def SetUp(self):
        twinplanes = (-2, 3)
        self.F = growth.Flake(*twinplanes)

    def test_prob(self):
        # probability distribution test
        pass


if __name__ == '__main__':
  unittest.main()
