#!/usr/bin/env python
## encoding: utf-8
"""                       FLAKE GROWTH SIMULATION
                          ~~~~~~~~~~~~~~~~~~~~~~~

                          Guilherme Stein © 2015
                          University of Würzburg
                          <guilherme.stein@physik.uni-wuerzburg.de>
"""
from __future__ import print_function, division, generators
import itertools as it
from random import choice
from helper import Grid
from mayavi import mlab as m


class Flake:
  """ Contains higher level methods for manipulating the flake.
  """
  def __init__(self, size=20, twins=(9, 11)):
    self.size = size
    self.twins = twins
    self.surface = []
    self.grid = Grid(size, twins)
    self.make_seed(radius=2)

  def get(self, idx):
    return self.grid.get(idx)

  def set(self, idx, **values):
    return self.grid.set(idx, **values)

  def clear(self):
    """ Erase all sites in the flake.
    """
    for site in self.permutator():
      self.grid.delete(site)


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


  def make_seed(self, radius=1):
    """ Populates the lattice points around the middle of the grid.

    This is used as initial flake seed.
    """
    self.clear()
    mid = self.size // 2
    for x in self.permutator(range(mid-radius, mid+radius+1)):
      self.grid.set(x, type='atom')
    self.create_surface()


  def occupied(self):
    """ Returns all sites in the grid that aren't empty.
    """
    return (site for site in self.permutator() if self.get(site))


  def chk(self, idx):
    for n in idx:
      if n not in range(self.size):
        return False
    return True


# # ########
#   grow  #
# #########
  def dim(self):
    D = []
    r = zip(*self.occupied())
    for i in r:
      delta = abs(max(i) - min(i))
      D.append(delta)
    return D


  def grow(self, rounds=1):
    for r in range(rounds):
      chosen = choice(self.surface)
      self.grid.set(chosen)
      self.change_surface(chosen)


  def create_surface(self):
    valid = (s for s in self.occupied() if self.get(s)['type'] == 'atom')
    self.surface = []
    for site in valid:
      self.surface_extender(site)


  def change_surface(self, idx):
    self.surface.remove(idx)
    self.surface_extender(idx)


  def surface_extender(self, idx):
    """ Extends the surface by the real neighbours of atom `idx`.
    """
    self.surface.extend([nb for nb, diffs in self.real_neighbours(idx) if not
                         self.grid.get(nb) and diffs < 2.3])


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
    absolutes = (tuple(sum(y) for y in x) for x in pairs)
    valids = (each for each in absolutes if self.chk(each))
    return valids


  def real_neighbours(self, idx):
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
    # print(choice_vec)
    all_indexed = list(self.abs_neighbours(idx))
    # print(all_indexed)
    all_coordinated = [self.grid.coord(site) for site in all_indexed]
    # print("AC:", all_coordinated)
    all_diffs = [choice_vec.dist(site) for site in all_coordinated]
    # print("AD:", all_diffs)
    all_associated = zip(all_indexed, all_diffs)
    aa = list(all_associated)
    return aa


# # ########
#   PLOT  #
# #########
  def plot(self, mayavi=True):
    """ Plot method of the flake.

    `scatter` expects three lists of xs, ys, zs, therefore the whole zip and
    unpack action. We only plot points that 'are' something.
    This is checked at creation of `valid` and then translated into vectors in
    `points`.
    If any argument is passed that evaluates `True` the surface will also be
    plotted.
    """
    _all = list(self.occupied()) + (not mayavi) * self.surface
    points = list(zip(*(self.grid.coord(site) for site in _all)))

    if mayavi:
      m.points3d(*points)
      m.show()
    else:
      import matplotlib.pyplot as plt
      from mpl_toolkits.mplot3d import Axes3D
      if False:
        Axes3D
      surface = self.surface

      def color(idx, t=None):
        if t == 'atom':
          color_list.append((1, 0, 0))
        elif t == 'surface':
          color_list.append((0.1, 0.1, 0.1))
        else:
          color_list.append((0.5, 0.5, 0.5))

      color_list = []
      for each in self.occupied():
        color(each, 'atom')
      for each in surface:
        color(each, 'surface')

      fig = plt.figure()
      ax = fig.add_subplot(111, projection='3d')
      ax.scatter(*points, s=1000, c=color_list)
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
  def avg(lst):
    added = sum(x for x in lst)
    return added/len(lst)
  runs = range(1)
  f = Flake(size=100, twins=(45, 50, 53, 66))
  print("Size :", f.size, "twins :", f.twins)
  dims = []
  steps = 10000
  for i in runs:
    f.make_seed()
    f.grow(steps)
    x = f.dim()
    # print(x)
    dims.append(x)
  ag = [sum(y)/len(runs) for y in zip(*dims)]
  print(ag)
  with open('../output/random_growth_01.dat', 'a') as fi:
    msg = ("\n-----------------------\n{s}::series\t"
    "@{st}::atoms\nAVERAGE (X Y Z):\t{av}").format(s=len(runs), st=steps, av=ag)
    fi.write(msg)
  # f.plot(surface=True)
  return dims, f


if __name__ == '__main__':
  main()
