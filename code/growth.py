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
import arrow
from os import path, makedirs

Q = False


class Flake:
  """ Contains higher level methods for manipulating the flake.
  """
  def __init__(self, size=20, twins=(), seed_size=1, height=10):
    self.size = size
    self.twins = twins
    self.height = height
    self.tag = ''
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
      self.set_surface(idx)
    else:
      self.grid.set(idx, value)

  def clear(self):
    for site in self.permutator():
      self.grid.delete(site)
    self.create_entire_surface()


# # ###########
#   GLOBALS  #
# ############
  def permutator(self, rng=None):
    """ Returns the 3D Cartesian product of `seed`.

    Without arguments it will return all sites in flake.
    """
    if not rng:
      return (i for i in it.product(range(self.size), repeat=3) if i[2] <
              self.height)
    else:
      return it.product(rng, repeat=3)


  def make_seed(self, radius=1):
    """ Populates the lattice points around the middle of the grid.

    This is used as initial flake seed.
    """
    mid = (self.size//2, self.size//2, self.height // 2)
    env = list(self.permutator(range(-radius, radius+1)))
    pairs = (zip(mid, others) for others in env)
    total = (tuple(sum(y) for y in x) for x in pairs)
    for each in total:
      self.set(each)
    self.create_entire_surface()


  def chk(self, idx):
    if idx[-1] not in range(self.height):
      return False
    for n in idx[:-1]:
      if n not in range(self.size):
        return False
    return True


# # #################
#   NEIGHBORHOOD  #
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
        qprint("Really NOTHING?! Found.", quiet=Q)

  def set_surface(self, site):
    occupied_neighbours = self.real_neighbours(site, void=False)
    slot = len(occupied_neighbours)
    try:
      self.surface[slot].append(site)
    except IndexError as e:
      qprint(e, "\nSite: ", site, "Too many empty neighbours: ", slot, quiet=Q)


  def create_entire_surface(self):
    self.surface = [[] for _ in range(12)]
    for atom in self.atoms:
      realz = self.real_neighbours(atom, void=True)
      for nb in realz:
        self.set_surface(nb)


# # ########
#   PLOT  #
# #########
  def plot(self, save=False, tag=''):
    """ The `color()` function appends values for colors to the `color_list`.
    This list must have the same length as `whole` to plot correctly.

    """
    surface_chain = list(it.chain.from_iterable(self.surface))

    whole = self.atoms + surface_chain
    x, y, z = list(zip(*(self.grid.coord(site) for site in whole)))

    def color(idx, t=None):
      if t == 'atom':
        color_list.append(0.8)
      elif t == 'surface':
        color_list.append(0.25)
      else:
        color_list.append(0.)

    color_list = []
    for each in self.atoms:
      color(each, 'atom')
    for each in surface_chain:
      color(each, 'surface')

    m.points3d(x, y, z, color_list, colormap="spectral", scale_factor=1.0,
                vmin=0, vmax=1.1)
    if save:
      m.options.offscreen = True
      if tag:
        tag = str(tag) + '_'
      save_dir = self.daily_output()
      _time = self.date[1].rsplit('.')[0].replace(':', '-')
      fname = path.join(save_dir, 'Flake@' + _time + '_S' + str(self.size) +
                        '_T' + str(self.twins) + tag + '.png')
      m.savefig(fname, size=(1024, 768))
    else:
      m.show()

  def daily_output(self):
    self.date = arrow.now().isoformat().rsplit('T')
    today = self.date[0]
    output_dir = path.join('../output/', today)
    output_dir = path.join(output_dir, self.tag)
    try:
      makedirs(output_dir)
    except OSError:
      pass
    return output_dir


# ----------
# -  main  -
# ----------
def main():
  f = Flake(size=31, twins=(), seed_size=0)
  counter = 1
  f.tag = 's3'
  for round in range(100):
    print("Flake now contains {} atoms.".format(counter))
    f.plot(save=True, tag=counter)
    counter += 50
    f.grow(20)


if __name__ == '__main__':
  main()
