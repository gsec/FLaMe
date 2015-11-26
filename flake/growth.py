#!/usr/bin/env python3
# encoding: utf-8
"""                       FLAKE GROWTH SIMULATION
                          ~~~~~~~~~~~~~~~~~~~~~~~

                          Guilherme Stein © 2015
                          University of Würzburg
                          <guilherme.stein@physik.uni-wuerzburg.de>
"""
from math import sqrt
import itertools as it
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from helper import Vector, GridError


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

# # ################
#   FUNDAMENTALS  #
# #################

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
    try:
      atom = self.raw_grid[i][j][k]
    except IndexError as e:
      print("An error occured while acessing a wrong index.", e)
      return None

    if not value:
      return atom[0]
    elif isinstance(value, (int, float)):
      self.raw_grid[i][j][k] = [True, value]
    elif value == 'get':
      return atom
    elif value == 'val':
      return atom[1]
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
    try:
      _perms = self.layer_permutations[k]
    except IndexError as e:
      print(e)
      return Vector(0, 0, 0)
    prototype = Vector(2*i + (j+k) % 2,
                       sqrt(3)*(j + _perms * 1/3),
                       k*2*sqrt(6)/3)
    return prototype

# # ###########
#   GLOBALS  #
# ############

  def whole_crystal(self):
    status = (self.grid(*idx, value='get') for idx in self.permutator())
    coords = (self.coord(*idx) for idx in self.permutator())
    return zip(status, coords)

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

  def create_surface(self):
    for atom in self.permutator():
      if self.grid(*atom, value='val') == 10:
        for nb in self.real_neighbours(*atom):
          if not self.grid(*nb):
            self.grid(*nb, value=2)

  def make_seed(self):
    i = self.size // 2
    for x in self.permutator((i-1, i, i+1)):
      self.grid(*x, value=10)

# # #################
#   NEIGHBOURHOOD  #
# ##################

  def real_neighbours(self, i, j, k, nn_switch='NN'):
    """ Creates next neighbours based on the distances.

    * build all (1, -1, 0) permutations
    * choose a site (ensure later that every possibility is covered...)
    * get coordinates of all surrounding sites
    * get vector differences of these sites with the chosen site
    * filter those, which are larger than two (two diameters equivalent to next
      neighbor)
    * zip them together

    """
    # TODO: not enough (<12) neighbours in some twinplane constellations

    choice = self.coord(i, j, k)
    all_indexed = self.abs_neighbours(i, j, k)
    all_coordinated = (self.coord(*site) for site in all_indexed)
    all_diffs = (choice.dist(site) for site in all_coordinated)
    all_associated = zip(all_indexed, all_diffs)
    if nn_switch == 'NN':
      next_relatives = [nb for nb, dist in all_associated if dist <= 2.1]
      # print('Atom {idx} with {vec} has {num} next neighbours.'.format(
          # idx=(i, j, k), vec=choice, num=len(next_relatives)))
      return next_relatives
    elif nn_switch == 'diff':
      next_diffs = [dist for nb, dist in all_associated if dist <= 2.1]
      return next_diffs

  def abs_neighbours(self, i, j, k):
    """ Return a list of next neighbours of `i, j, k`

    Here we should introduces a `idx_combinations` var and move (1, -1, 0)
    there.
    """
    relative_neighbours = self.permutator((1, -1, 0))
    relative_neighbours.remove((0, 0, 0))
    pairs = (zip((i, j, k), nn) for nn in relative_neighbours)
    return (tuple(sum(y) for y in x) for x in pairs)

# # ########
#   PLOT  #
# #########

  def plot(self, color=None):
    """ Plot method of the flake.

    `scatter` expects three lists of xs, ys, zs, therefore zip and unpack
    action. We only plot points that 'are' something.
    """
    whole = list(self.whole_crystal())
    status, pts = zip(*[(s, c) for s, c in whole if s[0]])
    if not color:
      color = []
      for each in status:
        if each[1] <= 1:
          color.append((1, 1, 1, 0))
        elif 1 < each[1] <= 2:
          color.append((0.5, 0.5, 0, 0.5))
        elif each[1] > 2:
          color.append((1, 0, 0, 1))
    points = list(zip(*pts))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(*points, s=1000, c=color)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    _rng = [0, 2 * self.size]
    ax.auto_scale_xyz(_rng, _rng, _rng)
    plt.show()


# ----------
# -  main  -
# ----------
def main():
  f = Flake(size=5, twins=(2, ))
  f.make_seed()
  f.create_surface()
  f.plot()


if __name__ == '__main__':
  main()
  Axes3D
