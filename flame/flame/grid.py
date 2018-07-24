"""
Utilities needed by the Flake Simulation.
Low level functions for `growth`
"""
from __future__ import print_function, division, generators
from math import sqrt
from itertools import product
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Grid(object):
    """ Basic Grid object.

    Here we keep track of the twin planes and according shifts.
    The `coords` method translates positions from index space into real space for this
    special grid configuration.
    """
    def __init__(self, twins):
        self.twins = twins
        self.twin_layers, self.upcounter, self.upsign = self.twin_gen()

    def twin_gen(self):
        """ Create a representation of the layer permutation .

        Return `list`, `counter`, `sign`

        Mapping an ABC layer to: A -> 0, B -> 1, C -> 2
        And flipping the order at each twin plane: ABCABCAB... -> AB'C'BACBA...
        where 'C' represents a twin plane. The `counter` keeps track of the permutation,
        while `sign` determines the permutation direction. Every twin plane inverts the
        permutation order, therefore `sign` is flipped at each twin plane we encounter.
        The generated layers are returned, along with the last `counter` and `sign` value
        to continue building the grid indefinitely.
        """
        if not self.twins:
            return None, None, None
        twin_layers = []
        sign = 1
        counter = (min(self.twins) - 1) % 3         # compensate for first addition
        twin_range = range(min(self.twins), max(self.twins) + 1)
        for layer in twin_range:
            counter += sign
            twin_layers.append(counter % 3)
            if layer in self.twins:
                sign *= -1
        return (twin_layers, counter % 3, sign)

    def shift(self, layer):
        """ Performs the permutation shift for twin planes TP = self.twin_gen().

        LOW:    Index is smaller than lowest twin plane, or there isn't a TP at all.
                Down to -iNf the order is as without twin planes. Just permutate happily
                mod 3.

        TP:     Index is in TP range, get shift directly.
                In the range of twin planes this is the list created by `twin_gen()`.

        HIGH:   Index is bigger than highest twinplane, shift index by `counter` as
                offset and continue in correct order given through `sign` up to +iNf.
        """
        if not self.twins or layer < min(self.twins):
            return layer % 3
        elif min(self.twins) <= layer <= max(self.twins):
            return self.twin_layers[layer - min(self.twins)]  # shifted layer index
        elif layer > max(self.twins):
            return (self.upcounter + self.upsign * (layer - max(self.twins))) % 3

    def coord(self, idx):
        """ Return Cartesian coordinates vector of a given lattice point.

        (i, j, k) are the indices and (a, b, c) are lattice base vectors. Crystal sites
        then are i*a + j*b + k*c.  `shift` is the permutation (0, 1 or 2) of the layer
        displacement according to the fcc-stacking and the twin plane configuration.
        """
        i, j, k = idx
        shift = self.shift(k)
        prototype = Vector(2*i + (j + shift) % 2,
                           sqrt(3)*(j + shift/3),
                           k*2*sqrt(6)/3)
        return prototype


class Seed():
    def __init__(self):
        self.seeds = {
                'point': set(
                    ((0, 0, 0),)
                ),
                'cube': set(
                    product((-1, 0, 1), repeat=3)
                ),
                'bigcube': set(
                    product((-2, -1, 0, 1, 2), repeat=3)
                ),
                'sphere': set(
                    ((0, -1, 1), (-1, 0, 1), (0, 0, 1),
                     (-1, -1, 0), (-1, 0, 0), (-1, 1, 0),
                     (0, 0, 0), (0, -1, 0), (0, 1, 0), (1, 0, 0),
                     (-1, -1, -1), (0, 0, -1), (0, -1, -1))
                ),
                'plane': set(
                    ((1, 0, 0), (1, 2, 0), (-2, 1, 0), (-2, 0, 0), (1, -1, 0),
                     (0, 1, 0), (-2, 2, 0), (-1, 0, 0), (-2, -2, 0), (0, -1, 0),
                     (1, 1, 0), (1, -2, 0), (0, -2, 0), (0, 2, 0), (2, 0, 0),
                     (-1, -2, 0), (-1, 1, 0), (-1, 2, 0), (2, -1, 0), (-1, -1, 0),
                     (2, 2, 0), (0, 0, 0), (2, 1, 0), (-2, -1, 0), (2, -2, 0)))
            }


    def seed_gen(self, shape=None):
        """ Create the first atoms to initialize the surface creation.

        `seeds` are sets of tuples containing atom indices. The return value is a
        tuple with the first element being the number of seed atoms and the second
        the set of tuples. If the requested shape is not found in the `SEEDS`
        dictionary it defaults to `point`.
        """
        try:
            seed = self.seeds[shape]
        except KeyError:
            if shape:
                logger.warn("Requested shape: {} not found! Defaulting to "
                            "`point`.".format(shape))
            seed = self.seeds['point']
        return len(seed), seed.copy()


class Vector(object):
    """ Self defined Vector object.

    Supports:
        * addition/subtraction with other vectors
        * addition/ subtraction with floats and integers (component wise)
        * equality if all components are equal
        * dist() method with other Vector as argument
        * length through the `abs()` method
    """
    def __init__(self, x, y, z):
        self.x, self.y, self.z = float(x), float(y), float(z)

    def __repr__(self):
        x, y, z = (round(r, 2) for r in (self.x, self.y, self.z))
        return 'Vector:({}, {}, {})'.format(x, y, z)

    def __iter__(self):
        """ Iteration over a Vector yields its components. """
        for comp in (self.x, self.y, self.z):
            yield comp

    def __eq__(self, other):
        """ Equality check is done by comparing each component. """
        eq_list = [self.__dict__[cmp] == other.__dict__[cmp] for cmp in ('x', 'y', 'z')]
        return all(eq_list)

    def __add__(self, other):
        """ Addition of vectors by component. """
        new_x = self.x + other.x
        new_y = self.y + other.y
        new_z = self.z + other.z
        return Vector(new_x, new_y, new_z)

    def __sub__(self, other):
        """ Subtraction of vectors by component. """
        new_x = self.x - other.x
        new_y = self.y - other.y
        new_z = self.z - other.z
        return Vector(new_x, new_y, new_z)

    def __mul__(self, other):
        """ Scaling of Vectors. """
        assert isinstance(other, (int, float)), "Vector scaling only with scalars!"
        new_x = other*self.x
        new_y = other*self.y
        new_z = other*self.z
        return Vector(new_x, new_y, new_z)

    def __rmul__(self, other):
        """ From both sides. """
        return self.__mul__(other)

    def __truediv__(self, other):
        """ Divides each component by a scalar. """
        assert isinstance(other, (int, float)), "Vector scaling only with scalars!"
        new_x = self.x/other
        new_y = self.y/other
        new_z = self.z/other
        return Vector(new_x, new_y, new_z)

    def dist(self, other):
        """ Return distance between the vectors. """
        delta = self - other
        return abs(delta)

    def __abs__(self):
        """ Return distance from origin. """
        return sqrt(self.x**2 + self.y**2 + self.z**2)
