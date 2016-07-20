"""
Utilities needed by the Flake Simulation.
"""
from __future__ import print_function, division, generators
from math import sqrt, pi


class Grid(object):
    """ Methods for getting the correct coordinates.

    Takes twin planes and the fcc-structure of the crystal into account.
    """
    def __init__(self, twins):
        self.twins = twins
        self.twin_layers, self.upcounter, self.upsign = self.twin_gen()


    def twin_gen(self):
        """ Create a representation of the layer permutation .

        Return `list`, `counter`, `sign`
        Mapping an ABC layer to: A -> 0, B -> 1, C -> 2
        And flipping the order at each twin plane: ABCABCAB... -> AB'C'BACBA...

        The `counter` keeps track of the permutation, while `sign` determines the
        permutation direction. `sign` is flipped at each twin plane we encounter.
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


    def shift(self, idx):
        """ Performs the permutation shift for twin planes TP=self.twin_gen().

        LOW:    Index is smaller than lowest twin plane, or there isn't a TP at all.
                Down to -iNf the order is as without twin planes. Just permutate happily
                mod 3.

        TP:     Index is in TP range, get shift directly.
                In the range of twin planes this is the list created by `twin_gen()`.

        HIGH:   Index is bigger than highest twinplane, shift index by `counter` as
                offset and continue in correct order given through `sign` up to +iNf.
        """
        if not self.twins or idx < min(self.twins):
            return idx % 3
        elif min(self.twins) <= idx <= max(self.twins):
            return self.twin_layers[idx - min(self.twins)]  # shifted layer index
        elif idx > max(self.twins):
            return (self.upcounter + self.upsign * (idx - max(self.twins))) % 3


    def coord(self, idx):
        """ Return Cartesian coordinates vector of a given lattice point.

        (i, j, k) are the indices and (a, b, c) are lattice base vectors. Crystal sites
        then are i*a + j*b + k*c.  `_perms` is the permutation (0, 1 or 2) of the layer
        displacement according to the fcc-stacking and the twin plane configuration. Every
        twin plane inverts the permutation order. The handling of the permutation is
        outsourced to `self.shift`.
        """
        i, j, k = idx
        shift = self.shift(k)
        prototype = Vector(2*i + (j + shift) % 2,
                           sqrt(3)*(j + shift/3),
                           k*2*sqrt(6)/3)
        return prototype


class Vector(object):
    """ Self defined Vector object.

    Supports:
        * addition/subtraction with other vectors
        * addition/ subtraction with floats and integers (for each component)
        * equality if components are equal
        * dist() method with other Vector as argument
        * length through the `abs()` method
    """
    def __init__(self, x, y, z):
        self.x, self.y, self.z = float(x), float(y), float(z)

    def __repr__(self):
        x, y, z = (round(r, 2) for r in (self.x, self.y, self.z))
        return 'Vector:({}, {}, {})'.format(x, y, z)


    def __iter__(self):
        """ Iteration over a Vector yields its components.  """
        for comp in (self.x, self.y, self.z):
            yield comp


    def __eq__(self, other):
        """ Equality check is done by comparing each component.  """
        eq_list = [self.__dict__[comp] ==
                             other.__dict__[comp] for comp in ('x', 'y', 'z')]
        return all(eq_list)


    def __add__(self, other):
        """ Addition of vectors by component.  """
        new_x = self.x + other.x
        new_y = self.y + other.y
        new_z = self.z + other.z
        return Vector(new_x, new_y, new_z)


    def __sub__(self, other):
        """ Subtraction of vectors by component.    """
        new_x = self.x - other.x
        new_y = self.y - other.y
        new_z = self.z - other.z
        return Vector(new_x, new_y, new_z)


    def __mul__(self, other):
        """ Scaling of Vectors.  """
        new_x = other*self.x
        new_y = other*self.y
        new_z = other*self.z
        return Vector(new_x, new_y, new_z)


    def __rmul__(self, other):
        return self.__mul__(other)


    def __div__(self, other):
        new_x = self.x/other
        new_y = self.y/other
        new_z = self.z/other
        return Vector(new_x, new_y, new_z)


    def dist(self, other):
        """ Return distance between the vectors.    """
        delta = self - other
        return abs(delta)


    def __abs__(self):
        """ Return distance from origin.
        """
        return sqrt(self.x**2 + self.y**2 + self.z**2)


class AtomsIO(object):
    """ Creates the atom object with two slots.
    """
    __slots__ = ('element', 'location')

    def __init__(self, element, location):
        self.element = element
        self.location = location
