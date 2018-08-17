## encoding: utf-8

import unittest
from os.path import isfile
from flame.growth import Flake


class TestFlakeBasics(unittest.TestCase):
    """ Flake sanity checks.

    Only for static flake tests.
    """
    def setUp(self):
        self.tF = Flake(seed='point')

    def test_seed_init(self):
        surface_one = set([(-1, 0, 1), (0, 0, 1), (0, -1, 1), (0, -1, 0),
                           (-1, 0, 0), (0, 1, 0), (-1, 1, 0), (1, 0, 0),
                           (-1, -1, 0), (0, -1, -1), (0, 0, -1), (-1, -1, -1)])
        self.assertEqual(self.tF.atoms, set([(0, 0, 0)]))
        self.assertEqual(self.tF.surface[1], surface_one)
        self.assertEqual(len(self.tF.integrated_surface()), len(self.tF.maxNB))

        self.assertEqual(self.tF.grid.twins, ())
        self.assertEqual(self.tF.grid.twin_layers, None)
        self.assertIsInstance(self.tF.__repr__(), str)


    def test_neighbours(self):
        """ Check if the chosen atom has the proper neighbours in index space.

        Ensure this translates into correct atoms considering the real space distance and
        make sure all neighbours return our seed as non-void neighborhood.
        """
        atom = (2, 5, 8)
        all_neighbours = [(2, 6, 9), (3, 4, 7), (3, 5, 9), (3, 4, 8), (1, 6, 7),
                          (3, 6, 7), (2, 4, 7), (2, 5, 7), (1, 5, 9), (1, 5, 8),
                          (3, 5, 8), (1, 5, 7), (2, 5, 9), (2, 4, 9), (2, 4, 8),
                          (3, 6, 9), (1, 6, 8), (3, 6, 8), (1, 6, 9), (1, 4, 9),
                          (3, 5, 7), (1, 4, 8), (3, 4, 9), (2, 6, 8), (2, 6, 7),
                          (1, 4, 7)]
        self.assertEqual(set(all_neighbours), set(self.tF.abs_neighbours(atom)))

        point_seed = (0, 0, 0)
        next_neighbours = [(0, -1, -1), (0, 0, -1),  (-1, 0, 1),  (-1, 0, 0),
                           (1, 0, 0),   (0, 0, 1),   (0, -1, 1),  (0, -1, 0),
                           (-1, 1, 0),  (-1, -1, 0), (0, 1, 0),   (-1, -1, -1)]

        self.assertEqual(set(next_neighbours),
                         set(self.tF.real_neighbours(point_seed, void=True)))

        for each in next_neighbours:
            self.assertEqual([point_seed], self.tF.real_neighbours(each, void=False))


    def test_twin_plane_creation(self):
        twinplanes = (4, 5, 6)
        f = Flake(*twinplanes)
        self.assertEqual(f.grid.twins, twinplanes)

        layers = [1, 0, 1]
        self.assertEqual(f.grid.twin_layers, layers)


class TestStaticFlake(unittest.TestCase):
    def setUp(self):
        twinplanes = (-2, 3)
        self.tF = Flake(*twinplanes)

    def test_prob_grow(self):
        self.tF.grow(mode='prob')

    def test_temp_dist(self):
        with self.assertRaises(IndexError):
            self.tF.temperature_dist(slot=-2, temp=10)

        with self.assertRaises(ValueError):
            self.tF.temperature_dist(slot=2, temp=-10)

        with self.assertRaises(ValueError):
            self.tF.temperature_dist(slot=2, temp=10e4)

        self.assertEqual(self.tF.temperature_dist(0), 0)

    def test_weights(self):
        simple = Flake(seed='point')
        initp = [0 for i in simple.maxNB]
        initp[1] = 12
        self.assertEqual(initp, simple.weights())


class TestFlakeGrowth(unittest.TestCase):
    def setUp(self):
        self.twins = (0, 1)
        self.seed = 'point'
        self.rounds = 37

    def test_base_modes(self):
        prbF = Flake(seed=self.seed, temp=0)
        detF = Flake(seed=self.seed)
        rndF = Flake(seed=self.seed)

        prbF.grow(self.rounds, mode='prob')
        detF.grow(self.rounds, mode='det')
        rndF.grow(self.rounds, mode='rand')

        for flake in (prbF, detF, rndF):
            self.assertEqual(flake.iter, self.rounds + 1)

    def test_carve(self):
        tF = Flake()
        tF.grow(self.rounds)
        tF.carve()

        self.assertGreater(tF.iter, len(tF.atoms))

    def test_cap_limit(self):
        detF = Flake()
        with self.assertRaises(StopIteration):
            detF.grow(self.rounds, mode='det', cap=10)

    def test_export_function(self):
        tF = Flake()
        tF.grow(self.rounds)
        fname = tF.export_coordinates('testexporter010')
        self.assertTrue(isfile(fname.name))

    def test_colorize(self):
        tF = Flake()
        self.assertIsInstance(tF.colorize(), dict)

    def test_plotting(self):
        """ Only testing the returned colors, not mayavi plotting...
        """
        tF = Flake()
        tF.grow(self.rounds)
        plot_data = list(tF.plot(ret=True))
        data_length = len(tF.atoms) + len(tF.integrated_surface())

        self.assertEqual(len(plot_data), 4)
        self.assertTrue(all(len(col) == data_length for col in plot_data))

    def test_bubbles(self):
        tF = Flake()
        while not tF.surface[0]:
            tF.grow(500, mode='rand')
