"""
Utilities needed by the Flake Simulation.
"""
from __future__ import print_function, division, generators
from math import sqrt


class Grid(object):
  def __init__(self, twins):
    """ Methods for getting the correct coordinates.

    Takes twin planes and the fcc-structure of the crystal into account.
    """
    self.twins = twins
    self.twin_layers, self.upcounter, self.upsign = self.twin_gen()


  def twin_gen(self):
    """ Create a z-list representing the permutation of the layer.

    Mapping an ABC layer to: A -> 0, B -> 1, C -> 2
    And flipping the order at each twin plane: ABCABCAB... -> AB'C'BACBA...
    """
    if not self.twins:
      return None, None, None
    twin_layers = []
    sign = 1
    counter = min(self.twins) % 3
    span = range(min(self.twins), max(self.twins) + 1)
    for layer in span:
      twin_layers.append(counter % 3)
      if layer in self.twins:
        sign *= -1
      counter += sign
    counter -= sign             # undo last addition for variable export
    return (twin_layers, counter % 3, sign)


  def shift(self, idx):
    """ Performs the twinplanes permutation shift.

    LOW:  Index is smaller than lowest twinplane, or there isn't a TP at all,
          just permutate happily mod 3
    TPs:  Index is in TP range, get shift from twin_gen() generated list
    HIGH: Index is bigger than highest twinplane, shift index by highest TP
          and continue in correct direction through sign
    """
    if not self.twins or idx < min(self.twins):
      return idx % 3
    elif min(self.twins) <= idx <= max(self.twins):
      return self.twin_layers[idx - min(self.twins)]  # shifted layer index
    elif idx > max(self.twins):
      return (self.upcounter + self.upsign * (idx - max(self.twins))) % 3


  def coord(self, (i, j, k)):
    """ Return Cartesian coordinates vector of a given lattice point.

    (i, j, k) are the indices and (a, b, c) are lattice base vectors. Crystal
    sites then are i*a + j*b + k*c.  `_perms` is the permutation (0, 1 or 2) of
    the layer displacement according to the fcc-stacking and the twin plane
    configuration. Every twin plane inverts the permutation order.
    """
    prototype = Vector((2*i + (j+k) % 2,
                      sqrt(3)*(j + self.shift(k) * 1/3),
                      k*2*sqrt(6)/3))
    return prototype


def qprint(*args, **kwargs):
  """ Override print function to suppress output.
  """
  try:
    if not kwargs.pop('quiet'):
      print(*args, **kwargs)
  except KeyError:
    pass


class Vector(object):
  """ Self defined Vector object.

  Supports:
    * addition/subtraction with other vectors
    * addition/ subtraction with floats and integers (for each component)
    * equality if components are equal
    * dist() method with other Vector as argument
    * length through the `abs()` method
  """
  def __init__(self, comp):
    self.comp = comp
    self.x, self.y, self.z = self.comp

  def __repr__(self):
    x, y, z = (round(r, 2) for r in self.comp)
    return 'Vector:({}, {}, {})'.format(x, y, z)


  def __iter__(self):
    """ Iteration over a Vector yields its components.
    """
    for comp in self.comp:
      yield comp


  def __eq__(self, other):
    """ Equality check is done by comparing each component.
    """
    eq_list = [self.__dict__[comp] ==
               other.__dict__[comp] for comp in ('x', 'y', 'z')]
    return all(eq_list)


  def __add__(self, other):
    """ Addition of vectors by component.

    Single numbers also can be added by applying it to each component
    separately.
    """
    try:
      new_x = self.x + other.x
      new_y = self.y + other.y
      new_z = self.z + other.z
    except:
      new_x = self.x + other
      new_y = self.y + other
      new_z = self.z + other
    return Vector((new_x, new_y, new_z))


  def __sub__(self, other):
    """ Subtraction of vectors by component.

    Single numbers also can be subtracted by applying it to each component
    separately.
    """
    try:
      new_x = self.x - other.x
      new_y = self.y - other.y
      new_z = self.z - other.z
    except:
      new_x = self.x - other
      new_y = self.y - other
      new_z = self.z - other
    return Vector((new_x, new_y, new_z))


  def dist(self, other):
    """ Require other Vector() object, return distance between the vectors.
    """
    if type(self) != type(other):
      msg = "Self:\t{s}\t{st}\nOther:\t{o}\t{ot}".format(
        s=self, st=type(self), o=other, ot=type(other))
      raise TypeError("Argument for `dist` must be another Vector." + msg)
    delta = self - other
    return abs(delta)


  def __abs__(self):
    """ Return distance from origin.
    """
    return sqrt(self.x**2 + self.y**2 + self.z**2)


class AtomsExport(object):
  """ Creates the atom object with two slots.
  """
  __slots__ = ('element', 'location')

  def __init__(self, element, location):
      self.element  = element
      self.location = location
