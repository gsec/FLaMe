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
from random import randrange
from helper import Grid


class Flake:
  """ Contains higher level methods for manipulating the flake.
  """

  def __init__(self, size=7, twins=()):
    self.size = size
    self.twins = twins
    self.grid = Grid(size, twins)


# # ###########
#   GLOBALS  #
# ############
  def permutator(self, seed=None):
    """ Creates all possible permutations of length three of all given objects
    in `seed`.

    Without arguments it will return all sites in flake.
    """
    if not seed:
      seed = range(self.size)
    types = it.combinations_with_replacement(seed, 3)
    perms = []
    for i in types:
      t = set(it.permutations(i))
      while t:
        perms.append(t.pop())
    return perms


  def clear(self):
    """ Erase all sites in the flake.
    """
    for site in self.permutator():
      self.grid.delete(site)


  def make_seed(self):
    """ Populates the lattice points around the middle of the grid.

    This is used as initial flake seed.
    """
    mid = self.size // 2
    for x in self.permutator((mid-1, mid, mid+1)):
      self.grid.set(x, type='atom', domain='seed')


  def create_surface(self):
    occupied = (s for s in self.permutator() if self.grid.get(s))
    valid = (s for s in occupied if self.grid.get(s)['type'] == 'atom')
    self.surface = []
    for site in valid:
      self.surface.extend(
        [nb for nb, diffs in self.real_neighbours(site) if not self.grid.get(nb) and diffs < 2.3])

    # atoms = []
    # for at in self.permutator():
      # try:
        # if self.grid.get(at)['type'] == 'atom':
          # atoms.append(at)
      # except TypeError:
        # pass
    # valid = (s for s in occupied if self.grid.get(s)['type'] == 'atom')
    # neigbours = (self.real_neighbours(s) for s in valid)
    # nbs = []
    # for nbtuple in neigbours:
      # val_idx = [(s, d) for s, d in nbtuple if self.size not in s]
      # val_type = [(s, d) for s, d in val_idx if not self.grid.get(s)['type'] == 'atom']
      # nbs.extend(val_type)
    # [self.grid.set(nb, type='surface') for nb, diffs in nbs if diffs < 3]
    # valid_neighbours = ((s, d) for s, d in neigbours if self.size not in s)
    # return neigbours
    # return occupied


# # ########
#   grow  #
# #########
  def grow(self, rounds=1):
    for r in range(rounds):
      self.create_surface()
      ra = randrange(len(self.surface))
      choice = self.surface[ra]
      self.grid.set(choice)
      # valid = [x for x in self.permutator() if self.grid.get(x)]
      # surf = [x for x in valid if self.grid.get(x)['type'] == 'surface']
      # return surf


# # #################
#   NEIGHBOURHOOD  #
# ##################
  def abs_neighbours(self, idx):
    """ Return a list of next neighbours of `i, j, k`

    Here we should introduces a `idx_combinations` var and move (1, -1, 0)
    there.
    """
    relative_neighbours = self.permutator((1, -1, 0))
    relative_neighbours.remove((0, 0, 0))
    pairs = (zip(idx, nn) for nn in relative_neighbours)
    absolutes = [tuple(sum(y) for y in x) for x in pairs]
    for v in absolutes:
      if not v:
        absolutes.remove(v)
        continue
      for c in v:
        if c not in range(self.size):
          absolutes.remove(v)
          break
    return absolutes


  def real_neighbours(self, idx, nn_switch='NN'):
    """ Creates next neighbours based on the distances.

    * choose a site (ensure later that every possibility is covered...)
    * build all (1, -1, 0) permutations
    * get coordinates of all surrounding sites
    * get vector differences between the atom and its neighbours
    * filter those, which are larger than two (two diameters equivalent to next
    neighbor)
      * zip them together and return zipped list as (indices, distance) pairs.
    """
    # TODO: Check for correct nn in different twin plane constellations.

    choice_vec = self.grid.coord(idx)
    all_indexed = self.abs_neighbours(idx)
    all_coordinated = (self.grid.coord(site) for site in all_indexed)
    all_diffs = (choice_vec.dist(site) for site in all_coordinated)
    all_associated = zip(all_indexed, all_diffs)
    return list(all_associated)


# # ########
#   PLOT  #
# #########
  def plot(self, color=None):
    """ Plot method of the flake.

    `scatter` expects three lists of xs, ys, zs, therefore the whole zip and
    unpack action. We only plot points that 'are' something.
    This is checked at creation of `valid` and then translated into vectors in
    `points`.
    """
    valid = [idx for idx in self.permutator() if self.grid.get(idx)]
    points = list(zip(*(self.grid.coord(site) for site in valid)))

    if not color:
      color = []
      for site, each in ((site, self.grid.get(site)['type']) for site in valid):
        if each == 'atom':
          color.append((1, 0, 0, 1))
        elif site in self.surface:
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
  f = Flake(size=8, twins=(2, 4))
  f.make_seed()
  f.create_surface()
  for times in range(54):
    f.create_surface()
    surf = f.grow()
    print(surf)
    luc = surf[randrange(0, len(surf))]
    f.grid.set(luc, type='atom')
    print(luc)
  f.plot()


if __name__ == '__main__':
  main()
  Axes3D
