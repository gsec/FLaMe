from math import pi
from collections import deque
from random import choice, random
import itertools as it
import logging

from flame.grid import Grid, Seed
from flame.settings import blender_helper, DIFF_CAP, AtomsIO

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Flake():
    """ Generates a Flake object with a FCC lattice and twin planes.

    The Flake take an arbitrary number of integers as `*args`, the twin planes.
    They are converted into a twin-plane tuple and proper layer stacking is
    calculated accordingly.

    Flake bootstrap
        The Flake is initialized with our (possibly empty) list of twin planes. The
        desired seed is evaluated and added to the atoms attribute. All settings have
        sane defaults for a typical growth. At the end we create the surface of
        possibilities around the newly generated atoms and weigh them according to the
        surrounding neighbours.

    Attr:
        maxNB (range): object we iterate over for neighbours or slots
            goes from 0 to 11 neighbours. This is because we can have zero neighbours
            if we fill a hole, but we can never have all 12 neighbours empty since we do
            not allow for single atoms detached from the Flake.

    Args:
        trail (int): trail length
            How many atoms will be marked as atoms trail

        seed (str): seed name
            from one of the following possibilities:[point, sphere, cube, bigcube, plane]

        temp (float): Artificial temperatures
            Accepted range [0 .. 1000].
            Lower temperatures lead to cleaner crystals in probabilty growth mode.
    """
    def __init__(self, *twins, **kwargs):
        self.twins = twins
        self.maxNB = range(12)

        self.seed_shape = kwargs.get('seed', 'sphere')
        self.iter, self.atoms = Seed().seed_gen(self.seed_shape)
        self.trail_length = kwargs.get('trail', 20)
        self.trail = deque(maxlen=self.trail_length)
        self.temp = kwargs.get('temp', 100)

        self.grid = Grid(twins)
        self._create_entire_surface()


    def __repr__(self):
        """ Representation of the Flake. [`twins`][`iterations`][`seed`][`temp`]
        """
        return ("fLakE :: twiNpLaNes:{}  Iterations:[{}]   Seed:[{}]  Temperature:"
                "[{}K]".format(self.twins, self.iter, self.seed_shape, self.temp))


####################
#     SURFACE     #
####################
    def _set_surface(self, site):
        """ Adds an empty site to the corresponding surface slot.
        """
        occupied_neighbours = self.real_neighbours(site, void=False)
        slot = len(occupied_neighbours)
        self.surface[slot].add(site)


    def _create_entire_surface(self):
        """ Generate the list for the surface sites based on occupied sites.

        Create a list of 12 sets (each corresponding to the  possibilities of next
        neighbours). Iterates through every atom and populates each adjacent empty site
        to the surface list.
        """
        self.surface = [set() for _ in self.maxNB]
        for atom in self.atoms:
            adjacent_voids = self.real_neighbours(atom, void=True)
            for nb in adjacent_voids:
                self._set_surface(nb)


    def integrated_surface(self):
        """ Return a set with all surface sites.
        """
        integrated_surface = set()
        for index in self.maxNB:
            integrated_surface.update(self.surface[index])
        return integrated_surface


    def sites(self):
        """ Return a list mapping binding slots to their quantity in surface.

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
        skin = set()
        srfc = self.integrated_surface()

        for spot in srfc:
            occupied_surface = self.real_neighbours(spot, void=False)
            skin.update(occupied_surface)

        self.atoms = skin


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


#########################
#     NEIGHBOURHOOD     #
#########################
    def abs_neighbours(self, atom):
        """ Return the list of next neighbours in index space of atom `idx`.
        """
        absNB = set(it.product(range(-1, 2), repeat=3))         # create NN indices
        absNB.remove((0, 0, 0))
        pairs = (zip(atom, nn) for nn in absNB)
        return [tuple(sum(y) for y in x) for x in pairs]


    def real_neighbours(self, atom, void=False):
        """ Creates next neighbours based on the distances.

        * build all neighbouring (1, -1, 0) index permutations of `atom`
        * get coordinates of all surrounding sites
        * get vector differences between the atom and its neighbours
        * filter those, which are larger than DIFF_CAP
        * return `void` or `occupied` neighbours

        TODO:  choose meaningful atomic radius
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
    def put_atom(self, at, slot):
        """ * remove atom from its surface slot
            * append to atoms list
            * check which neighbours are free, with those:
                * check every slot if it contains the neighbour
                * remove it from there
                * add it to next higher slot (cause now it has one more neighbour)
        """
        self.surface[slot].remove(at)
        self.atoms.add(at)

        if len(self.trail) >= self.trail.maxlen:
            self.trail.pop()
        self.trail.appendleft(at)       # prepends new atom to list of latest additions

        empty_neighbours = self.real_neighbours(at, void=True)
        for each in empty_neighbours:
            for e_slot, lst in enumerate(self.surface):
                if each in lst:
                    self.surface[e_slot].remove(each)
                    try:
                        self.surface[e_slot + 1].add(each)
                    except IndexError:
                        """ When this occours, we the surface is completely surrounded
                            by atoms, hence creating a bubble. We move those to surface
                            zero, which is otherwise unused.
                        """
                        self.surface[0].add(each)
                        logger.info("Bubble created...oO >> Site: {}".format(each))
                    break
            else:
                self.surface[1].add(each)        # create new surface entry for new ones
        self.iter += 1


    def grow(self, rounds=1, mode='prob', **kwargs):
        """ Here we manage which growth and other parameters as cap and iterations.

        The choices we get from each growth are then passed to the `put_atom` method,
        which then changes the atoms and surface of our flake.
        """
        modes = {'prob': self.prob_grow, 'rand': self.rand_grow, 'det': self.det_grow}

        for step in range(rounds):
            self.put_atom(*modes[mode](**kwargs))


    def prob_grow(self, rounds=1):
        """ The probabilistic growth mode, the most relevant in this simulation.

        The `temperature_dist` method returns the weights for each slot depending on
        slotnumber, temperature and number of atoms occupying that slot.

        All those weights are then added to our probability stack, where we randomly
        choose a point on this stack `stack_pointer`. To get the choice on which slot, we
        'consume' the weights-stack from the bottom upwards, until the 'stack-pointer' is
        reached. By breaking at this point we have the added probabilities as height on
        the stack and chosen our slot where the 'stack-pointer' lies.
        """
        # Here we save time by executing the function only once
        weights = self.weights()
        stack_pointer = random()*sum(weights)

        for slot, w in enumerate(weights):
            stack_pointer -= w
            if stack_pointer < 0:
                break
        chosen = choice(tuple(self.surface[slot]))

        return chosen, slot


    def weights(self):
        """ Return a list of probability weights for each slot in surface.

        Each element of the surface is weighted with the number of bindings to the
        power of a constant.
        """
        weights = []
        for slot in self.maxNB:
            p = self.temperature_dist(slot)
            weights.append(p)
        return weights


    def temperature_dist(self, slot, temp=None):
        """ Create the probabilty distribution for the surface.

        Depends on the slot (the more bindings,  the more probable it is to attach)
        scaled with an artificial temperature.
        """
        SCALE = 1/30            # scaling paramter
        UPPER_T = 1000          # temperature/numerical limits
        LOWER_T = 0

        if temp is None:
            temp = self.temp
        if not LOWER_T <= temp <= UPPER_T:
            raise ValueError("Temperature parameter must be between {} and {}.".format(
                LOWER_T, UPPER_T))
        if slot not in self.maxNB:
            raise IndexError("Slot must be in the range of next neighbours.")

        func = len(self.surface[slot]) * slot**(SCALE * (UPPER_T - temp))
        return func


    def det_grow(self, cap=0):
        """ Walk down the slots until the highest occupied slot is reached, then choose
        randomly an atom from there.
        """
        for slot in range(11, cap, -1):
            if self.surface[slot]:
                chosen = choice(tuple(self.surface[slot]))
                return chosen, slot
        else:
            raise StopIteration("Caplimit of <{}> for minimun free bindings reached.\n"
                                "Aborting further growth.".format(cap))


    def rand_grow(self):
        """ Choose uniformly between any valid surface site.
        """
        chosen = choice(list(it.chain.from_iterable(self.surface)))
        slot = next(i for i, x in enumerate(self.surface) if chosen in x)
        return chosen, slot


##################
#     OUTPUT     #
##################
    def export_coordinates(self, name):
        """ Exports the **xyz**-coordinates of the Flake atoms.
        Adapted from Atomic Blender, can be imported with xyz_io_mesh.
        Text format with a header and four columns: [ELEMENT, X, Y, Z]
        """
        geostr = ''.join(str((k, v)) + '\n' for k, v in self.geometry().items())
        header, xyzfile, infofile = blender_helper(name)

        attributes = header.format(twin=self.twins, temp=self.temp,
                                   shape=self.seed_shape, geo=geostr)

        output = ['{}\n'.format(len(self.atoms))]
        raw_atoms = (AtomsIO('Au', tuple(self.grid.coord(at))) for at in self.atoms)

        for atom in raw_atoms:
            string = '{:3s}{:15.5f}{:15.5f}{:15.5f}'.format(
                atom.element,
                atom.location[0],
                atom.location[1],
                atom.location[2])
            output.append(string)

        with open(xyzfile, 'w') as handler:
            handler.write('\n'.join(output))

        with open(infofile, 'w') as handler:
            handler.write(attributes)

        return infofile


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


    def plot(self, ret=False):
        """ Transforms the `colors` list to 4 lists of x, y, z, c.

        This list, the `clist` is either returned or displayed in mayavi.
        """
        coords = self.grid.coord
        colors = self.colorize()
        clist = zip(*[(coords(idx).x, coords(idx).y, coords(idx).z, c) for
                      (idx, c) in colors.items()])
        if ret:
            logger.info("Returning (x, y, z, colors) columns")
            return clist
        else:       # pragma: no cover
            try:
                from mayavi import mlab
                mlab.clf()
                mlab.points3d(*clist, colormap='gist_ncar',
                              scale_factor=0.1, vmin=0, vmax=15)
                mlab.show()
            except ImportError as e:
                logger.warn("Could not load {}. Is mayavi installed properly?\n"
                            "You may want to try flake.plot(ret=True) to return "
                            "the xyz coordinates.".format(e))
