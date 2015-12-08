"""
Utilities needed by the Flake Simulation.
"""
from __future__ import print_function, division, generators
from math import sqrt


class Vector(list):
  """ Self defined Vector object.

  Supports:
    * addition/subtraction with other vectors
    * addition/ subtraction with floats and integers (for each component)
    * length through the `abs()` method.
  """
  def __init__(self, x, y, z):
    self.x = x
    self.y = y
    self.z = z


  def __repr__(self):
    x, y, z = (round(r, 2) for r in (self.x, self.y, self.z))
    return 'Vector:({}, {}, {})'.format(x, y, z)


  def __iter__(self):
    """ Iteration over a Vector yields it's components.
    """
    for comp in (self.x, self.y, self.z):
      yield comp


  def __eq__(self, other):
    """ Equality check is done by comparing each component.
    """
    eq_list = [self.__dict__[comp] ==
               other.__dict__[comp] for comp in ('x', 'y', 'z')]
    return all(eq_list)


  def __add__(self, other):
    try:
      new_x = self.x + other.x
      new_y = self.y + other.y
      new_z = self.z + other.z
    except:
      new_x = self.x + other
      new_y = self.y + other
      new_z = self.z + other
    return Vector(new_x, new_y, new_z)


  def __sub__(self, other):
    try:
      new_x = self.x - other.x
      new_y = self.y - other.y
      new_z = self.z - other.z
    except:
      new_x = self.x - other
      new_y = self.y - other
      new_z = self.z - other
    return Vector(new_x, new_y, new_z)


  def dist(self, other):
    """ Require: Vector object. Return: distance between the vectors.
    """
    if not isinstance(other, self.__class__):
      raise TypeError("Argument for `dist` must be another Vector.")
    delta = self - other
    return abs(delta)


  def __abs__(self):
    return sqrt(self.x**2 + self.y**2 + self.z**2)


# # ###############
#   grid object  #
# ################

class Grid(list):
  def __init__(self, size, twins):
    """ Grid object containing the atom, accessed by the indices `i`, `j`, `k`.

    `size`: edge-length of the size**3 cube of lattice indices.
    `twins`: Iterable which yields numbers of twin plane layers.
    """
    self.size = size
    self.data = [[[None
      for _ in range(size)]
      for _ in range(size)]
      for _ in range(size)]
    self.layer_permutations = self.layer_gen(twins)


  def get(self, idx):
    i, j, k = idx
    try:
      return self.data[i][j][k]
    except IndexError:
      pass


  def set(self, idx, **value):
    i, j, k = idx
    if not self.data[i][j][k]:
      if not value:
        value = {'type': 'atom'}
      self.data[i][j][k] = value
    else:
      self.data[i][j][k].update(value)


  def delete(self, idx):
    i, j, k = idx
    self.data[i][j][k] = None


  def layer_gen(self, twins):
    """ Create a z-list representing the permutation of the layer.

    Mapping an ABC layer to: A -> 0, B -> 1, C -> 2
    And flipping the order at each twin plane: ABCABCAB... -> AB'C'BACBA...
    """
    L = []
    sign = 1
    counter = 0
    for layer in range(self.size):
      L.append(counter % 3)
      if layer in twins:
        sign = -1*sign
      counter += sign
    return L


  def coord(self, idx):
    """ Return Cartesian coordinates vector of a given lattice point.

    (i, j, k) are the indices and (a, b, c) are lattice base vectors. Crystal
    sites then are i*a + j*b + k*c.  `_perms` is the permutation (0, 1 or 2) of
    the layer displacement according to the fcc-stacking and the twin plane
    configuration. Every twin plane inverts the permutation order.
    """
    i, j, k = idx
    _perms = float(self.layer_permutations[k])
    prototype = Vector(2*i + (j+k) % 3,
                      sqrt(3)*(j + _perms * 1/3),
                      k*2*sqrt(6)/3)
    return prototype
