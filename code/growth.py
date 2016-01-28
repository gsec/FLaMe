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
import arrow
from random import choice
from helper import Grid, qprint
from mayavi import mlab as m
from os import path, makedirs

Q = False          # Set verbosity, Q is quiet


class AtomsExport(object):
  """ Creates the atom object with two slots.
  """
  __slots__ = ('element', 'location')

  def __init__(self, element, location):
      self.element  = element
      self.location = location


class Flake(object):
  """ Contains higher level methods for manipulating the flake.
  """
  def __init__(self, size=50, twins=None, seed_size=1, height=10):
    if twins is None:
      twins = (-1 + height//2, 1 + height//2)
    self.twins = twins
    self.size = size
    self.height = height
    self.seed_size = seed_size
    self.tag = ''
    self.iter = 0
    self.atoms = []
    self.surface = [[] for _ in range(12)]
    self.grid = Grid(size, twins, height)
    self.make_seed(radius=seed_size)

  def __repr__(self):
    return "fLakE ::\t sidelength[{}]\t twinplanes[{}]\t seed[{}]".format(
      self.size, self.twins, self.seed_size)

  def get(self, idx):
    """ Returns the integer at the index provided as tuple.
    """
    return self.grid.get(idx)

  def set(self, idx, value='atom'):
    """ Set value in `grid` and append to corresponding list attribute.
    """
    if value == 'atom' or 1:
      self.grid.set(idx, 1)
      self.atoms.append(idx)
    elif value == 'surface' or 2:
      self.grid.set(idx, 2)
      self.set_surface(idx)
    else:
      self.grid.set(idx, value)

  def clear(self):
    self.atoms = []
    self.surface = [[] for _ in range(12)]
    for site in self.permutator():
      self.grid.delete(site)


# # ###########
#   GLOBALS  #
# ############
  def permutator(self, rng=None):
    """ Returns the 3D Cartesian product of `seed`.

    Without arguments it will return all sites in flake.
    """
    ret = lambda R: (i for i in it.product(R, repeat=3) if self.chk(i))
    if not rng:
      return ret(range(self.size))
    else:
      return ret(rng)


  def make_seed(self, radius=1):
    """ Populates the lattice points around the middle of the grid.

    This is used as initial flake seed.
    """
    mid = (self.size//2, self.size//2, self.height // 2)
    env = list(self.permutator(range(-radius, radius+1)))
    pairs = (zip(mid, others) for others in env)
    total = (tuple(sum(y) for y in x) for x in pairs)
    for each in total:
      if self.chk(each):
        self.set(each)
    self.create_entire_surface()

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


  def chk(self, idx):
    """ Returns whether the index is in valid range.
    """
    if max(idx[0], idx[1]) >= self.size:
      return False
    elif idx[2] >= self.height:
      return False
    else:
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
    if not self.chk(atom):
      return []
    choice_vec = self.grid.coord(atom)
    indexed = list(self.abs_neighbours(atom))
    coordinates = (self.grid.coord(ai_site) for ai_site in indexed)
    diffs = (choice_vec.dist(ac_site) for ac_site in coordinates)
    associated = list(zip(indexed, diffs))
    nearest = [nb for nb, diff in associated if diff < DIFF_CAP]
    if void:
      return [nb for nb in nearest if not self.grid.get(nb)]
    else:
      return [nb for nb in nearest if self.grid.get(nb)]


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
          self.surface[slot].remove(chosen)
          self.set(chosen)
          his_neighbours = self.real_neighbours(chosen, void=True)
          for site in his_neighbours:
            realz = self.real_neighbours(site)
            bindings = len(realz)
            if site not in self.surface[bindings]:
              self.surface[bindings].append(site)
          self.iter += 1
          break
      else:
        qprint("Really NOTHING?! Found.", quiet=Q)


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

    m.clf()
    m.points3d(x, y, z, color_list, colormap="spectral",
                         scale_factor=1.0, vmin=0, vmax=1.1)
    if save:
      m.options.offscreen = True
      if tag:
        tag = str(tag) + '_'
      save_dir = self.daily_output()
      _time = self.date[1].rsplit('.')[0].replace(':', '-')
      fname = path.join(save_dir, 'Flake@' + _time + '_S' + str(self.size) +
                        '_T' + str(self.twins) + tag + '.png')
      m.savefig(fname, size=(1024, 768))
      m.close()
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

  def export(self, tag=''):
    """ Simplified export function adapted from 'io_mesh_xyz'.
    """
    raw_atoms = (('Au', tuple(self.grid.coord(at)))
                 for at in self.atoms if at)
    list_atoms = []
    counter = 0

    for each in raw_atoms:
      list_atoms.append(AtomsExport(*each))
      counter += 1

    if self.tag and not tag:
      tag = self.tag
    elif not tag and not self.tag:
      tag = 'flake'
    save_dir = self.daily_output()
    fname = "{tag}.{width}x{height}_TP-{tp}_it-{it}_{tag}.xyz".format(
      width=self.size, height=self.height, tp=self.twins, it=len(self.atoms),
      tag=tag)
    filepath_xyz = path.join(save_dir, fname)
    with open(filepath_xyz, "w") as xyz_file_p:
      xyz_file_p.write("%d\n" % counter)
      xyz_file_p.write("This is a XYZ file according to Atomic Blender - XYZ. "
                      "***WITH MODIFICATIONS! TAKE CARE AND READ THE CODE*** "
                      "For more details see: wiki.blender.org/index.php/"
                      "Extensions:2.6/Py/Scripts/Import-Export/XYZ\n")

      for i, atom in enumerate(list_atoms):
          string = "%3s%15.5f%15.5f%15.5f\n" % (
                                        atom.element,
                                        atom.location[0],
                                        atom.location[1],
                                        atom.location[2])
          xyz_file_p.write(string)


# ----------
# -  main  -
# ----------
def animate(tag, binning=20):
  f = Flake(size=71, seed_size=0)
  atomic_num = (2*f.seed_size + 1)**3
  f.tag = tag

  all_timings = '\n'
  # f.plot(save=True, tag='seed')
  for round in range(500):
    atomic_num += binning
    g_start = arrow.now()
    f.grow(binning)
    g_end = arrow.now()
    p_start = arrow.now()
    f.plot(save=True, tag=atomic_num)
    p_end  = arrow.now()
    p_delta = (p_end - p_start).total_seconds()
    g_delta = (g_end - g_start).total_seconds()
    timing_string = ("Count: {} atoms. Added {} atoms in {} "
                     "sec and rendered in {} sec\n"
                     ).format(atomic_num, binning, g_delta, p_delta)
    qprint(timing_string, quiet=Q)
    all_timings += timing_string
  fname = path.join(f.daily_output(), 'timings_' + f.tag + '.txt')
  with open(fname, 'a+') as tfile:
    tfile.write(f.__repr__() + all_timings + '\n')


def main():
  animate('lateNtrial', binning=20)
  # for x in range(10):
    # tag = 'pro002_single_' + str(x)
    # binn = (1 + x) * 50
    # animate(tag, binning=1)

if __name__ == '__main__':
  main()
