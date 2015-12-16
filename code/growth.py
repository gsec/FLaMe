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
from helper import Grid, qprint
from mayavi import mlab as m

Q = False


class Flake:
  """ Contains higher level methods for manipulating the flake.
  """
  def __init__(self, size=20, twins=(), seed_size=1, height=10):
    self.size = size
    self.twins = twins
    self.height = height
    self.atoms = []
    self.surface = [[] for _ in range(12)]
    self.grid = Grid(size, twins, height)
    self.make_seed(radius=seed_size)

  def get(self, idx):
    return self.grid.get(idx)

  def set(self, idx, value='atom'):
    if value == 'atom' or 1:
      self.grid.set(idx, 1)
      self.atoms.append(idx)
    elif value == 'surface' or 2:
      self.grid.set(idx, 2)
    else:
      self.grid.set(idx, value)

  def clear(self):
    for site in self.permutator():
      self.grid.delete(site)
    self.create_entire_surface()


# # ###########
#   GLOBALS  #
# ############
  def permutator(self, seed=None):
    """ Returns the 3D cartesian product of `seed`.

    Without arguments it will return all sites in flake.
    """
    if not seed:
      return (i for i in it.product(range(self.size), repeat=3) if i[2] <
              self.height)
    else:
      return it.product(seed, repeat=3)

  # def permutator(self, seed=None):
    # """ Creates all possible permutations of length three of all given objects
    # in `seed`.

    # Without arguments it will return all sites in flake.
    # """
    # if not seed:
      # seed = range(self.size)
    # types = it.combinations_with_replacement(seed, 3)
    # perms = []
    # for i in types:
      # t = set(it.permutations(i))
      # while t:
        # perms.append(t.pop())
    # return perms


  def make_seed(self, radius=1):
    """ Populates the lattice points around the middle of the grid.

    This is used as initial flake seed.
    """
    mid = self.height // 2
    for x in self.permutator(range(mid-radius, mid+radius+1)):
      self.set(x)
    self.create_entire_surface()


  def occupied(self):
    """ Returns all sites in the grid that aren't empty.
    """
    # TODO: substitute self.occupied() with self.atoms list
    return (site for site in self.permutator() if self.get(site))


  def chk(self, idx):
    for n in idx:
      if n not in range(self.size):
        return False
    return True


# # #################
#   NEIGHBOURHOOD  #
# ##################
  def abs_neighbours(self, idx):
    """ Return a list of next neighbours of `i, j, k`

    Here we should introduces a `idx_combinations` var and move (1, -1, 0)
    there.
    """
    relative_neighbours = list(self.permutator((1, -1, 0)))
    relative_neighbours.remove((0, 0, 0))
    pairs = (zip(idx, nn) for nn in relative_neighbours)
    absolutes = (tuple(sum(y) for y in x) for x in pairs)
    valids = (each for each in absolutes if self.chk(each))
    return valids


  def real_neighbours(self, atom, void=False):
    """ Creates next neighbours based on the distances.

    * choose a site (ensure later that every possibility is covered...)
    * build all (1, -1, 0) permutations
    * get coordinates of all surrounding sites
    * get vector differences between the atom and its neighbours
    * filter those, which are larger than DIFF_CAP
    * return `void` or `occupied` neighbours
    """
    DIFF_CAP = 2.3
    choice_vec = self.grid.coord(atom)
    all_indexed = list(self.abs_neighbours(atom))
    all_coordinated = [self.grid.coord(site) for site in all_indexed]
    all_diffs = [choice_vec.dist(site) for site in all_coordinated]
    all_associated = zip(all_indexed, all_diffs)
    aa = list(all_associated)
    only_next = [nb for nb, diffs in aa if diffs < DIFF_CAP]
    if void:
      return [nb for nb in only_next if not self.grid.get(nb)]
    else:
      return [nb for nb in only_next if self.grid.get(nb)]


# # ########
#   grow  #
# #########
  def grow(self, rounds=1):
    """ Manipulates the flake, adding atoms to its surface.
    """
    for r in range(rounds):
      for slot in range(11, 0, -1):
        if self.surface[slot]:
          chosen = choice(self.surface[slot])
          i = 0
          while chosen in self.surface[slot]:     # kinda hacked.... !!
            qprint(i, quiet=True)
            i += 1
            self.surface[slot].remove(chosen)
          self.set(chosen)
          his_neighbours = self.real_neighbours(chosen, void=True)
          for site in his_neighbours:
            realz = self.real_neighbours(site)
            bindings = len(realz)
            self.surface[bindings].append(site)
          break
      else:
        qprint("Really NOTHING?! found.", quiet=Q)

  def set_surface(self, site):
    occupied_neighbours = self.real_neighbours(site, void=False)
    slot = len(occupied_neighbours)
    try:
      self.surface[slot].append(site)
    except IndexError as e:
      qprint(e, "\nSite: ", site, "Too many empty neighbours: ", slot, quiet=Q)



  def create_entire_surface(self):
    atoms = (s for s in self.occupied() if self.get(s) == 1)
    self.surface = [[] for _ in range(12)]
    for atom in atoms:
      realz = self.real_neighbours(atom, void=True)
      for nb in realz:
        self.set_surface(nb)


# # ########
#   PLOT  #
# #########
  def plot(self):
    """ The `color()` function appends values for colors to the `color_list`.
    This list must have the same length as `whole` to plot correctly.

    """
    sfc = list(it.chain.from_iterable(self.surface))
    occ = list(self.occupied())
    whole = occ + sfc
    x, y, z = list(zip(*(self.grid.coord(site) for site in whole)))

    def color(idx, t=None):
      if t == 'atom':
        color_list.append(0.8)
      elif t == 'surface':
        color_list.append(0.25)
      else:
        color_list.append(0.)

    color_list = []
    for each in occ:
      color(each, 'atom')
    for each in sfc:
      color(each, 'surface')

    m.points3d(x, y, z, color_list, colormap="spectral", scale_factor=1.0,
                vmin=0, vmax=1.1)
    m.show()


# ----------
# -  main  -
# ----------
def main():
  f = Flake(size=31, twins=(), seed_size=0)
  f.grow(1000)
  f.plot()


if __name__ == '__main__':
  main()
