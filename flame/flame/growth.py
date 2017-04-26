#!/usr/bin/env python
## encoding: utf-8

from __future__ import print_function, division, generators
from math import pi
from arrow import now
from collections import deque
from os import path, makedirs
from random import choice, random
import itertools as it
import logging

from flame.grid import Grid
from flame.settings import (GROW_OUTPUT, PICKLE_EXT, DIFF_CAP,
                            AtomsIO, get_time, seed_gen)


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# This handles the `pickle` module for python 2 and 3
try:
    import cPickle as pickle
except ImportError:
    import pickle


class Flake(object):
    """ Generates a Flake object with a FCC lattice and twin planes.

    The Flake take an arbitrary number of integers as `*args`, the twin planes.
    They are converted into a twin-plane tuple and proper layer stacking is
    calculated accordingly.

    Following keyword arguments are possible:
        `seed`: str in {point, sphere, cube, bigcube, plane}
        `trail`: int for length of the marked atoms trail
        `temp`: float in {0 .. 273}, determines the probability through exponential.
    """
    DATE = get_time()


    @staticmethod
    def daily_output(name):
        """ Return output folder according to date.

        Create the folder if not already existent.
        """
        if len(name.split('/')) == 1:
            new_path = path.join(GROW_OUTPUT, Flake.DATE[0], name)
        else:
            new_path = path.join(GROW_OUTPUT, name)
        if not path.exists(path.dirname(new_path)):
            makedirs(path.dirname(new_path))
        return new_path


    @staticmethod
    def load(name):
        """ Return a flake instance from pickled file.
        """
        fname = Flake.daily_output(name) + PICKLE_EXT
        with open(fname, 'rb') as file_handler:
            logger.info("Loading Flake instance from file {} ...".format(fname))
            return pickle.load(file_handler)


    def save(self, name):
        """ Save the flake as pickled instance.

        For further analysis, continuation of growth with different parameters and simple
        pause/continue of growth.
        """
        fname = Flake.daily_output(name) + PICKLE_EXT
        with open(fname, 'wb') as file_handler:
            pickle.dump(self, file_handler, protocol=2)
        logger.info("Flake instance saved to disk: {}".format(fname))


    def __init__(self, *twins, **kwargs):
        """ Flake bootstrap.

        The Flake is initialized with our (possibly empty) list of twin planes. The
        desired seed is evaluated and added to the atoms attribute. All settings have
        sane defaults for a typical growth. At the end we create the surface of
        possibilities around the newly generated atoms and weigh them according to the
        surrounding neighbors.

        `maxNB` goes from 0 to 11 neighbors. This is because we can have zero neighbors
        if we fill a hole, but we can never have all 12 neighbors empty since we do not
        allow for single atoms detached from the Flake.
        """
        self.twins = twins
        self.maxNB = range(12)

        self.seed_shape = kwargs.get('seed', 'sphere')
        self.iter, self.atoms = seed_gen(self.seed_shape)
        self.trail_length = kwargs.get('trail', 20)
        self.trail = deque(maxlen=self.trail_length)
        self.temp = kwargs.get('temp', 50)

        self.grid = Grid(twins)
        self.create_entire_surface()


    def __repr__(self):
        """ Representation of the Flake. [`twins`][`iterations`][`seed`][`temp`]
        """
        return ("fLakE :: twiNpLaNes:{}  Iterations:[{}]   Seed:[{}]  Temperature:"
                "[{}K]".format(self.twins, self.iter, self.seed_shape, self.temp))


  ####################
  #     SURFACE     #
  ####################
    def set_surface(self, site):
        """ Adds an empty site to the corresponding surface slot.
        """
        occupied_neighbours = self.real_neighbours(site, void=False)
        slot = len(occupied_neighbours)
        self.surface[slot].add(site)


    def create_entire_surface(self):
        """ Generate the list for the surface sites based on occupied sites.

        Create a list of 12 sets (each corresponding to the  possibilities of next
        neighbors). Iterates through every atom and populates each adjacent empty site
        to the surface list.
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
        for index in self.maxNB:
            integrated_surface.update(self.surface[index])
        return integrated_surface


    def sites(self):
        """ Return a dictionary mapping binding slots to their quantity in surface.

        This allows to check quickly for the distribution of free binding sites on the
        surface.
        """
        return [len(x) for x in self.surface]


    def carve(self):
        """ Makes the flake hollow by removing atoms without adjacent surface.

        Optimizes the visualization of the Flake when plotted with mayavi. Growth can
        continue normally, since the inward empty space in the Flake is not part of the
        surface. This will fail if the `create_entire_surface()` method is used as it
        would also create an inner surface.
        """
        t_start = now()
        skin = set()
        srfc = self.integrated_surface()

        for spot in srfc:
            occupied_surface = self.real_neighbours(spot, void=False)
            skin.update(occupied_surface)
        diff = len(self.atoms) - len(skin)
        self.atoms = skin

        t_end = now()
        t_delta = (t_end - t_start).total_seconds()
        logger.info("Carved out {} atoms in {} sec.".format(diff, t_delta))


    def geometry(self):
        """ Calculate geometry information about the Flake.

        Here we generate all the data of the growth process we will later use. The
        dictionary returned here will be the data we collect in growth simulations.
        We can add other data with the dict.update() method, but existing fields should
        not be changed for compatibility reasons.
        `height`: Distance between furthest atoms in z direction (+1, taking the outside)
        `radius`: distance from center to far most atom
        `area`: pi*radius**2
        `aspect ratio`: comparison between the vertical and horizontal dimension
                                        (2*radius/height)
        """
        COORD = self.grid.coord
        POOL = self.atoms
        SITES = self.sites()

        mxz = max(POOL, key=lambda i: i[2])[2]
        mnz = min(POOL, key=lambda i: i[2])[2]
        max_outer = max(POOL, key=lambda i: i[0]**2 + i[1]**2)
        mean_binds = sum(weight*i for i, weight in
                         enumerate(num/sum(SITES) for num in SITES))

        attr = {
            'radius': COORD(max_outer).dist(COORD((0, 0, max_outer[2])))
        }
        attr.update({                # + 1 correct for border
            'height': COORD((0, 0, mxz + 1)).dist(COORD((0, 0, mnz)))
        })
        attr.update({
            'aspect_ratio': 2*attr['radius']/attr['height']
        })
        attr.update({
            'area': pi*attr['radius']**2
        })
        attr.update({
            'layers': mxz + 1 - mnz
        })
        attr.update({
            'iter': self.iter
        })
        attr.update({
            'bindings': mean_binds
        })


        for (k, i) in attr.items():
            setattr(self, k, i)
        return attr


  ########################
  #     NEIGHBORHOOD     #
  ########################
    def abs_neighbours(self, atom):
        """ Return the list of next neighbors in index space of atom `idx`.
        """
        absNB = set(it.product(range(-1, 2), repeat=3))         # create NN indices
        absNB.remove((0, 0, 0))
        pairs = (zip(atom, nn) for nn in absNB)
        return [tuple(sum(y) for y in x) for x in pairs]


    def real_neighbours(self, atom, void=False):
        """ Creates next neighbors based on the distances.

        * choose a site (ensure later that every possibility is covered...)
        * build all (1, -1, 0) permutations
        * get coordinates of all surrounding sites
        * get vector differences between the atom and its neighbors
        * filter those, which are larger than DIFF_CAP
        * return `void` or `occupied` neighbors
        """
        indexed = self.abs_neighbours(atom)
        coordinates = (self.grid.coord(ai_site) for ai_site in indexed)
        choice_vec = self.grid.coord(atom)
        diffs = (choice_vec.dist(ac_site) for ac_site in coordinates)
        associated = zip(indexed, diffs)

        nearest = (nb for nb, diff in associated if diff < DIFF_CAP)
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
            logger.debug("[{}]  :{}:  +{} rounds   CAP {}".format(
                mode.upper(), self.iter, rounds, cap))
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
            logger.debug("WEIGHTS: {}  Sum:{:e}".format(weights, sum(weights)))
            rnd = random()*sum(weights)
            for slot, w in enumerate(weights):
                logger.debug("RAND:{:e}  Slot:{}  Weight:{:e}".format(rnd, slot, w))
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
                logger.warn("Flake has no sites with {} or more free bindings. "
                                 "STOP at {}th growth step.\n".format(cap, step))
                ans = input("Continue with single growth? [y/n]  ")
                if ans in 'Yy':
                    return True
                elif ans in 'nN':
                    return False

        t_start = now()
        go()
        t_delta = (now() - t_start).total_seconds()
        return t_delta


    def put_atom(self, at, slot):
        """ * remove atom from its slot
            * append it to atoms list
            * check which neighbors are free, with those:
                * check every slot if it contains the neighbor
                * remove it from there
                * add it to next higher slot (cause now it has one more neighbor)
        """
        self.surface[slot].remove(at)
        self.atoms.add(at)
        if len(self.trail) >= self.trail.maxlen:
            self.trail.pop()
        self.trail.appendleft(at)          # prepends new atom to list of latest additions

        empty_neighbours = self.real_neighbours(at, void=True)
        for each in empty_neighbours:
            for e_slot, lst in enumerate(self.surface):
                if each in lst:
                    self.surface[e_slot].remove(each)
                    try:
                        self.surface[e_slot + 1].add(each)
                    except IndexError:
                        logger.debug("Filled a bubble...oO Site: {}".format(each))
                    break
            else:
                self.surface[1].add(each)        # create new surface entry for new ones
        self.iter += 1


  ##################
  #     OUTPUT     #
  ##################
    def export(self, name):
        """ Exports the **xyz**-coordinates of the Flake atoms.
        Adapted from Atomic Blender, can be imported with xyz_io_mesh.
        Text format with a header and four columns: [ELEMENT, X, Y, Z]
        """
        raw_atoms = (AtomsIO('Au', tuple(self.grid.coord(at))) for at in
                     self.atoms)

        with open(Flake.daily_output(name) + '.xyz', 'w') as xyz_file:
            xyz_file.write('{}\n{}\nThis is a XYZ file.\n'.format(
              len(self.atoms), str(self.geometry())))
            for atom in raw_atoms:
                string = '{:3s}{:15.5f}{:15.5f}{:15.5f}\n'.format(
                  atom.element,
                  atom.location[0],
                  atom.location[1],
                  atom.location[2])
                xyz_file.write(string)

    def colorize(self):
        """ Generate colors for visual representation of atoms and surface sites.

        Assign each atom and surface element a fix value according to their
        slot position, i.e. the amount of binding options. We have special
        values for the center atom and the trail, the last atoms added.
        Returns dict.
        """
        logger.info("Generating colors...")
        colors = dict((at, 15) for at in self.atoms)
        colors.update((at, 13) for at in self.trail)
        colors.update({(0, 0, 0): 14})      # mark the middle spot
        colors.update([(entry, slot) for slot, shelf in enumerate(self.surface)
                       for entry in shelf])
        logger.info("Color generation finished.")
        return colors


    def plot(self, save=False, pipeline=False):
        """ Transforms the `colors` list to 4 lists of x, y, z, c.

        This list, the `clist` is either returned or displayed in mayavi.
        Optionally a picture is saved to disk.
        """
        coords = self.grid.coord
        colors = self.colorize()
        clist = zip(*[(coords(idx).x, coords(idx).y, coords(idx).z, c) for
                      (idx, c) in colors.items()])
        if pipeline:
            # used for external mayavi rendering
            return clist

        try:
            from mayavi import mlab
        except ImportError as e:
            logger.warn("Module {} not available for Python 3. Flake can not be plotted."
                        "Returning (x, y, z, colors) columns".format(str(e).split()[-1]))
            return clist

        mlab.clf()
        mlab.points3d(*clist, colormap='gist_ncar', scale_factor=0.1, vmin=0, vmax=15)

        if not save:
            mlab.show()
        else:
            mlab.options.offscreen = True      # currently not working (bug in mayavi?)
            fname = 'Flake@' + Flake.DATE[1] + '_T' + str(self.twins) + '.png'
            fpath = Flake.daily_output(fname)
            mlab.savefig(fpath, size=(1024, 768))
            logger.info("Figure saved to {}".format(fpath))
            mlab.close()
