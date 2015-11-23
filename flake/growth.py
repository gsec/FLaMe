#!/usr/bin/env python3
# encoding: utf-8
"""                       FLAKE GROWTH SIMULATION
                          ~~~~~~~~~~~~~~~~~~~~~~~

                          Guilherme Stein         : 2015
                          University of Würzburg
                          <guilherme.stein@physik.uni-wuerzburg.de>
"""
from math import sqrt
import itertools as it
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from helper import Vector


# ---------------
# -  The Flake  -
# ---------------
class Flake:
  """ Class containing the atom and lattice informations about the flake."""

  def __init__(self, size=7, twins=tuple()):
    """
    Initialize `raw_grid` as size**3 array of [0, 0].
    """
    self.size = size
    self.twins = twins
    self.layer_permutations = self.layer_gen(*twins)
    self.raw_grid = [[[
                    [False, 0]
                    for _ in range(size)]
                    for _ in range(size)]
                    for _ in range(size)]

  def grid(self, i, j, k, value=None):
    """ Interface to access the `grid_list` containing information about each
    atom.

    Needs indices (i,j,k) and can access the slots `occupation` and `energy` of
    the atom on this site. The `raw_grid` has  a list of length two for each
    atom it represents at each site of the lattice. The first one is a boolean
    that signals if there is an atom at this site, the second one is a float or
    int representing the energy and therefore the associated probability of an
    atom attaching to this site. If called only with the indices it returns True
    or False, whether there is an atom present at the site. The three index
    arguments are mandatory. Further there is an optional *value* parameter
    which can be:
      'set'   : Sets the occupation entry to True.
      'get'   : Returns the occupation and energy for that site.
      'remove': Deletes atom from this site, setting `occupation` to False and
                `energy` to zero.
      #float  : Set `energy` entry to that value.
    """
    atom = self.raw_grid[i][j][k]

    if not value:
      return atom[0]
    elif isinstance(value, (int, float)):
      self.raw_grid[i][j][k] = [True, value]
    elif value == 'get':
      return atom
    elif value == 'set':
      self.raw_grid[i][j][k][0] = True
    elif value == 'remove':
      self.raw_grid[i][j][k] = [False, 0]
    else:
      raise GridError("Option `" + str(value) + "` not available.")

  def coord(self, i, j, k):
    """ Return Cartesian coordinates vector of a given lattice point.

    (i, j, k) are the indices and (a, b, c) are lattice base vectors. Crystal
    sites then are i*a + j*b + k*c.  `_perms` is the permutation (0, 1 or 2) of
    the layer displacement according to the fcc-stacking and the twin plane
    configuration. Every twin plane inverts the permutation order.
    """
    _perms = self.layer_permutations[k]
    prototype = Vector(2*i + j % 2 + k % 3,
                       sqrt(3)*(j + _perms * 1/3),
                       k*2*sqrt(6)/3)
    return prototype

  def layer_gen(self, *twins):
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

  def real_neighbours(self, i, j, k):
    """ Creates next neighbours based on the distances.

    * build all (1, -1, 0) permutations
    * choose a site (ensure later that every possibility is covered...)
    * get coordinates of all surrounding sites
    * get vector differences of these sites with the chosen site
    * filter those, which are larger than two (two diameters equivalent to next
      neighbor)
    * zip them together

    Performs about one order of magnitude slower than self.neighbours(), while
    containing it as a step.
    """
    # TODO: not enough (<12) neighbours in some constellations/permutations
    # CHECK THAT !!!

    choice = i, j, k
    all_relative = self.permutator((1, -1, 0))
    all_relative.remove((0, 0, 0))              # (0, 0, 0) is not a neighbor
    all_indexed = list(self.neighbours(*choice,
                                       relative_neighbours=all_relative))
    all_coordinated = [self.coord(*site) for site in all_indexed]
    # round for numerical issues
    all_diffs = [site.dist(self.coord(*choice)) for site in all_coordinated]
    all_associated = list(zip(all_relative, all_indexed, all_diffs))
    next_relatives = [each[0] for each in all_associated if each[2] <= 2.1]
    print('len of nn is: ', len(next_relatives))
    return next_relatives

  def neighbours(self, i, j, k, relative_neighbours=None):
    """ Return a list of next neighbours of `i, j, k` """
    site = (i, j, k)
    if not relative_neighbours:
      relative_neighbours = self.real_neighbours(*site)
    pairs = (zip(site, nn) for nn in relative_neighbours)
    return (tuple(sum(y) for y in x) for x in pairs)
    # return [Vector(*x) for x in ret]

  def set_neighbours(self, i, j, k, val=None):
    for n in self.neighbours(i, j, k):
      if not val:
        self.grid(*n, value='set')
      else:
        self.grid(*n, value=val)

  def create_surface(self):
    self.surface = []
    for atom in self.permutator():
      for nb in self.real_neighbours(*atom):
        if nb not in self.surface:
          self.surface.extend(nb)

  def permutator(self, seed=None):
    """ Creates all possible permutations of length three of all given objects
    in `seed`. """
    if not seed:
      seed = range(self.size)
    types = it.combinations_with_replacement(seed, 3)
    perms = []
    for i in types:
      t = set(it.permutations(i))
      while t:
        perms.append(t.pop())
    return perms

  def whole_crystal(self):
    return (self.coord(*idx) for idx in self.permutator())

  def plot(self, sites=None, color=['yellow']):
    """ Plot method of the flake. """
    if not sites:
      sites = self.whole_crystal()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    pts = list(zip(*sites))
    ax.scatter(*pts, s=1000, c=color)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    _rng = [0, 2 * self.size]
    ax.auto_scale_xyz(_rng, _rng, _rng)
    plt.show()


class GridError(Exception):
  pass


class InsufficientNeighbours(GridError):
  pass


class LayerError(GridError):
  pass


# ----------
# -  main  -
# ----------
def main(i=2, j=2, k=2):
  f = Flake(size=5, twins=(2, ))
  special_one = (i, j, k)
  coords = []
  cols = []
  f.grid(*special_one, value=5)
  f.set_neighbours(*special_one, val=3)
  print('\nRendering:')
  for idx in f.permutator(range(f.size)):
    _c = f.coord(*idx)
    coords.append(_c)
    cols.append(f.grid(*idx, value='get')[1] + 0.1 * idx[2])
  f.plot(coords, cols)


if __name__ == '__main__':
  main()
  Axes3D
