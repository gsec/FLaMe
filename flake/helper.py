"""
Utilities needed but not especially part of the Flake Project.
"""
from math import sqrt


class Vector:
  """ Self defined Vector object.

  Supports addition/subtraction with other vectors, addition/ subtraction with floats and
  integers (for each component) and the length through the `abs()` method.
  """
  def __init__(self, x, y, z):
    self.x = x
    self.y = y
    self.z = z

  def __repr__(self):
    return 'Vector:({}, {}, {})'.format(self.x, self.y, self.z)

  def __iter__(self):
    """ Iteration over a Vector yields it's components. """
    for comp in (self.x, self.y, self.z):
      yield comp

  def __eq__(self, other):
    eq_list = [self.__dict__[comp] == other.__dict__[comp] for comp in ('x', 'y', 'z')]
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
    if not isinstance(other, self.__class__):
      raise TypeError("Argument for `dist` must be another Vector.")
    delta = self - other
    return abs(delta)

  def __abs__(self):
    return sqrt(self.x**2 + self.y**2 + self.z**2)


class VectorException(Exception):
  pass
