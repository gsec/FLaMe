#!/usr/bin/env python3
# encoding: utf-8
"""                       FLAKE GROWTH SIMULATION
                          ~~~~~~~~~~~~~~~~~~~~~~~

                          Guilherme Stein © 2015
                          University of Würzburg
                          <guilherme.stein@physik.uni-wuerzburg.de>
"""
import itertools as it
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from helper import Grid


# ---------------
# -  The Flake  -
# ---------------
class Flake:
  """ Class containing the atom and lattice informations about the flake."""

  def __init__(self, size=7, twins=()):
    """
    # Initialize `raw_grid` as size**3 array of [0, 0].

    Here we will replace each spot on the grid with a NONE value.
    """
    self.size = size
    self.grid = Grid(size, twins)

# # ###########
#   GLOBALS  #
# ############

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
    occupied = (s for s in self.permutator() if self.grid.get(s))
    valid = (s for s in occupied if self.grid.get(s)['type'] == 'atom')
    for site in valid:
      [self.grid.set(nb, type='surface') for nb, diffs in
       self.real_neighbours(site) if not self.grid.get(nb)]

  def make_seed(self):
    i = self.size // 2
    for x in self.permutator((i-1, i, i+1)):
      self.grid.set(x, type='atom')

# # #################
#   NEIGHBOURHOOD  #
# ##################

  def real_neighbours(self, idx, nn_switch='NN'):
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

    choice_vec = self.grid.coord(idx)
    all_indexed = self.abs_neighbours(idx)
    all_coordinated = (self.grid.coord(site) for site in all_indexed)
    all_diffs = (choice_vec.dist(site) for site in all_coordinated)
    all_associated = zip(all_indexed, all_diffs)
    return list(all_associated)

  def abs_neighbours(self, idx):
    """ Return a list of next neighbours of `i, j, k`

    Here we should introduces a `idx_combinations` var and move (1, -1, 0)
    there.
    """
    relative_neighbours = self.permutator((1, -1, 0))
    relative_neighbours.remove((0, 0, 0))
    pairs = (zip(idx, nn) for nn in relative_neighbours)
    return (tuple(sum(y) for y in x) for x in pairs)

# # ########
#   PLOT  #
# #########

  def plot(self, color=None):
    """ Plot method of the flake.

    `scatter` expects three lists of xs, ys, zs, therefore zip and unpack
    action. We only plot points that 'are' something.
    """
    valid = [idx for idx in self.permutator() if self.grid.get(idx)]
    points = list(zip(*(self.grid.coord(site) for site in valid)))

    if not color:
      color = []
      for each in (self.grid.get(site)['type'] for site in valid):
        if each == 'atom':
          color.append((1, 0, 0, 1))
        elif each == 'surface':
          color.append((0.1, 0.1, 0.1, 0.3))
        else:
          color.append((0.5, 0.5, 0, 0.5))

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
