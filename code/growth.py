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
from random import choice, randrange
from helper import Grid, qprint
from mayavi import mlab as m


class Flake:
  """ Contains higher level methods for manipulating the flake.
  """
  def __init__(self, size=20, twins=()):
    self.size = size
    self.twins = twins
    self.grid = Grid(size, twins)
    self.make_seed(radius=1)
    self.create_surface()

  def get(self, idx):
    return self.grid.get(idx)

  def set(self, idx, **values):
    return self.grid.set(idx, **values)

  def clear(self):
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


  def real_neighbours(self, atom, only_void=False):
    """ Creates next neighbours based on the distances.

    * choose a site (ensure later that every possibility is covered...)
    * build all (1, -1, 0) permutations
    * get coordinates of all surrounding sites
    * get vector differences between the atom and its neighbours
    * filter those, which are larger than two (two diameters equivalent to next
    neighbor)
      * zip them together and return zipped list as (indices, distance) pairs.
    """
    DIFF_CAP = 2.3
    choice_vec = self.grid.coord(atom)
    all_indexed = list(self.abs_neighbours(atom))
    all_coordinated = [self.grid.coord(site) for site in all_indexed]
    all_diffs = [choice_vec.dist(site) for site in all_coordinated]
    all_associated = zip(all_indexed, all_diffs)
    aa = list(all_associated)
    only_next = [nb for nb, diffs in aa if  diffs < DIFF_CAP]
    if only_void:
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
          qprint("Accessed Slot #: ", slot, "\tSlot content: ",
                self.surface[slot], "\tAtomic choice: ", chosen)
          # tmp_list = self.surface[slot][:]
          # print(tmp_list)
          # tmp_list.remove(chosen)
          # print(tmp_list)
          # self.surface[slot] = tmp_list
          # print(self.surface[slot])
          while chosen in self.surface[slot]:
            self.surface[slot].remove(chosen)
            # raise Exception("NOT REMOVING CHOICE FROM SURFACE")
          self.set(chosen)
          his_neighbours = self.real_neighbours(chosen, only_void=True)
          for site in his_neighbours:
            realz = self.real_neighbours(site)
            bindings = len(realz)
            if bindings > 11:
              qprint("ALERT\n")
              qprint("site", site)
              qprint("bind-num", bindings)
              qprint("real neighbours of site:", realz)
              qprint("\n\n")
              self.surface[11].append(site)
              break
            self.surface[bindings].append(site)
          break
      else:
        qprint("Really NOTHING?! found.")


  def create_surface(self):
    atoms = (s for s in self.occupied() if self.get(s)['type'] == 'atom')
    self.surface = [[] for _ in range(12)]
    for atom in atoms:
      realz = self.real_neighbours(atom, only_void=True)
      for nb in realz:
        nb2nb = self.real_neighbours(nb)
        binds = len(nb2nb)
        if nb not in self.surface[binds]:
          try:
            self.surface[binds].append(nb)
          except IndexError as e:
            # self.surface[11].append(nb)
            qprint(e)
            qprint("Site: ", nb, "Has too many possible neighbours: ", binds)
            qprint("This must be an empty site")


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
    sfc = list(it.chain.from_iterable(self.surface))
    occ = list(self.occupied())

    def color(idx, t=None):
      if t == 'atom':
        color_list.append(randrange(98, 100) / 100)
      elif t == 'surface':
        color_list.append(randrange(25, 30) / 100)
      else:
        color_list.append(0.)

    color_list = []
    for each in self.occupied():
      color(each, 'atom')
    for each in sfc:
      color(each, 'surface')

    whole = occ + sfc
    qprint("CLIST", color_list)

    x, y, z = list(zip(*(self.grid.coord(site) for site in whole)))
    qprint("Xyz:", x)
    qprint('xlen', len(x), '\nclistlen', len(color_list))

    if mayavi:
      m.points3d(x, y, z, color_list, colormap="spectral", scale_factor=.8,
                 vmin=0, vmax=1.1)
      m.show()
      return occ, sfc
    else:
      import matplotlib.pyplot as plt
      from mpl_toolkits.mplot3d import Axes3D
      if False:
        Axes3D
      fig = plt.figure()
      ax = fig.add_subplot(111, projection='3d')
      ax.scatter(x, y, z, s=1000, c=color_list)
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
  f = Flake(size=11, twins=())
  f.plot()


if __name__ == '__main__':
  main()
