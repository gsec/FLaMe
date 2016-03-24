#!/usr/bin/env python2
## encoding: utf-8
"""                           FLaMe - a FLakeLAtticeMOdelEr

                              <guilherme.stein@physik.uni-wuerzburg.de>
"""
from __future__ import print_function, division, generators
import arrow
import cPickle
import itertools as it
import logging
from helper import *
from collections import deque
from os import path, mkdir, environ
from mayavi import mlab as m
from random import choice, random


class Flake(object):
    """ Generates a whole Flake object with a FCC lattice and twin planes.

    The Flake take an arbitrary number of integers as `*args`, the twin planes.
    They are converted into a twin-plane tuple and proper layer stacking is
    calculated accordingly.

    Following keyword arguments are possible:
        `seed`: str in {point, sphere, cube, bigcube, plane}
        `trail`: int for length of the marked atoms trail
        `temp`: float in {0 .. 273}, determines the probability through exponential.
    """

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    OUTPUT_DIR = path.join(environ['THESIS_PATH'], 'output')
    DATE = arrow.now().isoformat().rsplit('T')
    DIFF_CAP = 2.1


    @staticmethod
    def daily_output(fpath=None):
        """ Generator for the date output folder.

        Creates the folder if not already existent.
        """
        if not fpath:
            fpath = path.join(Flake.OUTPUT_DIR, 'TEMP')
        elif len(fpath.split('/')) == 1:
            fpath = path.join(Flake.OUTPUT_DIR, Flake.DATE[0], fpath)
        else:
            fpath = path.join(Flake.OUTPUT_DIR, fpath)
        if not path.exists(path.dirname(fpath)):
            mkdir(path.dirname(fpath))
        return fpath


    @staticmethod
    def load(*args):
        """ Return a flake instance from pickled file.
        """
        with open(Flake.daily_output(*args) + '.flm', 'rb') as file_handler:
            return cPickle.load(file_handler)

    def __init__(self, *twins, **kwargs):
        """ Flake bootstrap.

        * Create a seed
        * Add those to atoms list
        * Generate appropriate surface
        """
        self.twins = twins
        self.maxNB = range(12)                          # take all 12 NN as possibilities

        self.seed_shape = kwargs.get('seed', 'sphere')
        self.trail_length = kwargs.get('trail', 20)
        self.temp = kwargs.get('temp', 50)
        self.iter, self.atoms = self.seed(self.seed_shape)
        self.trail = deque(maxlen=self.trail_length)    # specify length of trail

        self.grid = Grid(twins)
        self.create_entire_surface()
        self.colors = self.color_init()


    def __repr__(self):
        """ Representation of the Flake. [`twins`][`iterations`][`seed`][`temp`]
        """
        return ("fLakE :: twiNpLaNes:{}  Iterations:[{}]   Seed:[{}]  Temperature:"
                "[{}K]".format(self.twins, self.iter, self.seed_shape, self.temp))


  ##########################
  #     INITIALIZATION     #
  ##########################
    def seed(self, shape='point'):
        """ Create the first atoms to initialize the surface creation.

        `seeds` are sets of tuples containing atom indices. The return value is a
        tuple with the first element being the number of seed atoms and the second
        the set of tuples.
        """
        if shape == 'point':
            seed = set(((0, 0, 0),))
        elif shape == 'cube':
            seed = set(it.product((-1, 0, 1), repeat=3))
        elif shape == 'bigcube':
            seed = set(it.product((-2, -1, 0, 1, 2), repeat=3))
        elif shape == 'sphere':
            seed = set(((0, -1, 1), (0, 0, 1), (-1, 0, 1), (-1, -1, -1), (0, -1, -1),
                        (0, 0, -1), (-1, -1, 0), (-1, 0, 0), (-1, 1, 0), (0, -1, 0),
                        (0, 1, 0), (1, 0, 0), (0, 0, 0)))
        elif shape == 'plane':
            seed = set(((1, 0, 0), (1, 2, 0), (-2, 1, 0), (-2, 0, 0), (1, -1, 0),
                        (0, 1, 0), (-2, 2, 0), (-1, 0, 0), (-2, -2, 0), (0, -1, 0),
                        (1, 1, 0), (1, -2, 0), (0, -2, 0), (0, 2, 0), (2, 0, 0),
                        (-1, -2, 0), (-1, 1, 0), (-1, 2, 0), (2, -1, 0), (-1, -1, 0),
                        (2, 2, 0), (0, 0, 0), (2, 1, 0), (-2, -1, 0), (2, -2, 0)))
        return len(seed), seed

    def color_init(self):
        """ Initialize the `color` attribute.

        Assign each atom a fixed value and surface elements according to their slot
        position. This is stored as dict in `self.colors`. The colors-attribute
        keeps track of atoms and surface values for visual representation.
        """
        color_dict = dict((at, 15) for at in self.atoms)
        color_dict.update([(entry, slot) for slot, shelf in enumerate(self.surface)
                                             for entry in shelf])
        return color_dict


  ####################
  #     SURFACES     #
  ####################
    def set_surface(self, site):
        """ Adds an empty site to the corresponding surface slot.
        """
        occupied_neighbours = self.real_neighbours(site, void=False)
        slot = len(occupied_neighbours)
        self.surface[slot].add(site)


    def create_entire_surface(self):
        """ Create a list of 12 sets (NN possibilities) and populates them.

        Iterates through every atom and populates each adjacent empty site to the
        surface list.
        """
        self.surface = [set() for _ in self.maxNB]
        for atom in self.atoms:
            adjacent_voids = self.real_neighbours(atom, void=True)
            for nb in adjacent_voids:
                self.set_surface(nb)


    def integrated_surface(self):
        """ Return a set with all surface sites.
        """
        integrated_surface = set()
        for i in self.maxNB:
            integrated_surface.update(self.surface[i])
        return integrated_surface


    def sites(self):
        """ Return a dictionary mapping binding slots to their quantity in surface.
        """
        return {i: len(x) for (i, x) in enumerate(self.surface)}


    def carve(self):
        """ Makes the flake hollow by removing atoms without adjacent surface.

        Affects self.atoms and self.colors
        """
        t_start = arrow.now()
        skin = set()
        srfc = self.integrated_surface()

        for spot in srfc:
            occupied_surface = self.real_neighbours(spot, void=False)
            skin.update(occupied_surface)
        diff = len(self.atoms) - len(skin)
        self.atoms = skin

        whole = srfc.union(skin)
        new_dict = {x: self.colors[x] for x in whole}
        self.colors = new_dict

        t_end = arrow.now()
        t_delta = (t_end - t_start).total_seconds()
        Flake.logger.info("Carved out {} atoms in {} sec.".format(diff, t_delta))


    def geometry(self):
        """ Calculate geometry information about the Flake.

        `height`: diff + 1 of furthest atoms in z-direction
        `radius`: distance from center to farmost atom
        `area`: pi*radius**2
        `aspect ratio`: comparison between the vertikal and horizontal dimension
                                        (2*radius/height)
        """
        COORD = self.grid.coord
        POOL = self.atoms

        mxz = max(POOL, key=lambda i: i[2])[2]
        mnz = min(POOL, key=lambda i: i[2])[2]
        mxr = max(POOL, key=lambda i: i[0]**2 + i[1]**2)

        items = {                                    # + 1 correct for border
            'height': COORD((0, 0, mxz + 1)).dist(COORD((0, 0, mnz))),
            'radius': COORD(mxr).dist(COORD((0, 0, mxr[2])))}
        items.update({
            'area': pi*items['radius']**2,
            'aspect_ratio': 2*items['radius']/items['height'],
            'layers': mxz + 1 - mnz})

        for (k, i) in items.iteritems():
            setattr(self, k, i)
        return items


  ########################
  #     NEIGHBORHOOD     #
  ########################
    def abs_neighbours(self, atom):
        """ Return the list of next neighbours in index space of atom `idx`.
        """
        absNB = set(it.product(range(-1, 2), repeat=3))         # create NN indices
        absNB.remove((0, 0, 0))
        pairs = (zip(atom, nn) for nn in absNB)
        return [tuple(sum(y) for y in x) for x in pairs]


    def real_neighbours(self, atom, void=False):
        """ Creates next neighbours based on the distances.

        * choose a site (ensure later that every possibility is covered...)
        * build all (1, -1, 0) permutations
        * get coordinates of all surrounding sites
        * get vector differences between the atom and its neighbours
        * filter those, which are larger than DIFF_CAP
        * return `void` or `occupied` neighbours
        """
        indexed = self.abs_neighbours(atom)
        coordinates = (self.grid.coord(ai_site) for ai_site in indexed)
        choice_vec = self.grid.coord(atom)
        diffs = (choice_vec.dist(ac_site) for ac_site in coordinates)
        associated = zip(indexed, diffs)

        nearest = (nb for nb, diff in associated if diff < Flake.DIFF_CAP)
        if void:
            return [nb for nb in nearest if nb not in self.atoms]
        else:
            return [nb for nb in nearest if nb in self.atoms]


  ##################
  #     GROWTH     #
  ##################
    def prob(self):
        """ Return a list of probability weights for each slot in surface.

        Each element of the surface is weighted with the number of bindings to the
        power of a constant.
        """
        weights = []
        func = lambda x: len(self.surface[x]) * x**(abs(self.temp))
        for slot in self.maxNB:
            p = func(slot)
            weights.append(p)
        return weights


    def grow(self, rounds=1, mode='prob', cap=1):
        """ Transform a surface site into an atom.

        `rounds` is the number of growth iterations. We have three different modes:
        'prob', 'det' and 'rand'.  If there are no more surface sites with equal or more
        than a `cap` number of bindings, growth will stop; then ask the user to continue.
        """
        def go():
            Flake.logger.info("GROWTH PROCEDURE -=-=-=- [{}]  {} rounds   "
                             "CAP {}".format(mode.upper(), rounds, cap))
            func_dict = {'prob': prob_grow, 'det': det_grow, 'rand': rand_grow}
            for step in range(int(rounds)):
                func = func_dict[mode]
                while not any(self.sites()[x] for x in range(cap, 12)):
                    if ask(step):
                        break
                    else:
                        return None
                self.put_atom(*func())

        def prob_grow():
            weights = self.prob()
            Flake.logger.debug("WEIGHTS: {}  Sum:{:e}".format(weights, sum(weights)))
            rnd = random()*sum(weights)
            for slot, w in enumerate(weights):
                Flake.logger.debug("RAND:{:e}  Slot:{}  Weight:{:e}".format(rnd, slot, w))
                rnd -= w
                if rnd < 0:
                    break
            chosen = choice(tuple(self.surface[slot]))
            return chosen, slot

        def det_grow():
            for slot in range(11, 0, -1):
                if self.surface[slot]:
                    chosen = choice(tuple(self.surface[slot]))
                    return chosen, slot

        def rand_grow():
            chosen = choice(list(it.chain.from_iterable(self.surface)))
            slot = next(i for i, x in enumerate(self.surface) if chosen in x)
            return chosen, slot

        def ask(step):
            while True:
                Flake.logger.warn("Flake has no sites with {} or more free bindings. "
                                 "STOPP at {}th growth step.\n".format(cap, step))
                ans = raw_input("Continue with single growth? [y/n]  ")
                if ans in 'Yy':
                    return True
                elif ans in 'nN':
                    return False

        t_start = arrow.now()
        go()
        t_delta = (arrow.now() - t_start).total_seconds()
        return t_delta


    def put_atom(self, at, slot):
        """ * remove atom from its slot
            * append it to atoms list
            * check which nb are free, with those:
                * check every slot if it contains the nb
                * remove it from there
                * add it to next higher slot (cause now he's got one nb more)
        """
        self.surface[slot].remove(at)
        self.atoms.add(at)
        self.colors.update(((at, 13),))
        if len(self.trail) >= self.trail.maxlen:
            old = self.trail.pop()
            self.colors.update(((old, 15),))
        self.trail.appendleft(at)          # prepend new atom to list of latest additions

        empty_neighbours = self.real_neighbours(at, void=True)
        for each in empty_neighbours:
            for e_slot, lst in enumerate(self.surface):
                if each in lst:
                    self.surface[e_slot].remove(each)
                    try:
                        self.surface[e_slot + 1].add(each)
                        self.colors.update(((each, e_slot + 1),))
                    except IndexError:
                        self.colors.pop(each)
                        Flake.logger.warn("Filled a bubble...oO")
                    break
            else:
                self.surface[1].add(each)        # create new surface entry for new ones
                self.colors.update(((each, 1),))
        self.iter += 1


  ##################
  #     OUTPUT     #
  ##################
    def save(self, *args):
        """ Saves the flake as pickled instance.
        """
        with open(Flake.daily_output(*args) + '.flm', 'wb') as file_handler:
            cPickle.dump(self, file_handler)

    def export(self, *args):
        """ Exports the **xyz**-coordinates of the Flake atoms.
        Adapted from Atomic Blender, can be imported with xyz_io_mesh.
        Text format with a header and four columns: [ELEMENT, X, Y, Z]
        """
        raw_atoms = (AtomsIO('Au', tuple(self.grid.coord(at))) for at in
                     self.atoms)

        with open(Flake.daily_output(*args) + '.xyz', 'w') as xyz_file:
            xyz_file.write('{}\n{}\nThis is a XYZ file.\n'.format(
              len(self.atoms), str(self.geometry())))
            for atom in raw_atoms:
                string = '{:3s}{:15.5f}{:15.5f}{:15.5f}\n'.format(
                  atom.element,
                  atom.location[0],
                  atom.location[1],
                  atom.location[2])
                xyz_file.write(string)


    def plot(self, save=False, pipeline=False):
        """ Transforms the `colors` list to 4 lists of x, y, z, c.

        This list, the `clist` is either returned or displayed in mayavi.
        Optionally a picture is saved to disk.
        """
        co = self.grid.coord
        clist = zip(*[(co(idx).x, co(idx).y, co(idx).z, c) for (idx, c) in
                                    self.colors.items()])
        if pipeline:
            return clist

        m.clf()
        m.points3d(*clist, colormap="gist_ncar", scale_factor=0.1, vmin=0, vmax=15)
        if save == 1:
            m.options.offscreen = False      # currently not working (bug in mayavi?)
            _time = Flake.DATE[1].rsplit('.')[0].replace(':', '-')
            fname = 'Flake@' + _time + '_T' + str(self.twins) + '.png'
            fpath = Flake.daily_output(fname)
            m.savefig(fpath, size=(1024, 768))
            m.close()
        else:
            m.show()
