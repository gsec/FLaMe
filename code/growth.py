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
from random import choice, random
from math import pi
from helper import *
from mayavi import mlab as m
from os import path, mkdir

Q = False          # Set verbosity, Q is quiet


class Flake(object):
  """ FLAKE creation class.
  """
  def __init__(self, *twins):
    self.iter = 0
    self.twins = twins
    self.span = range(12)                 # takes all 12 NN as possibilities
    self.atoms = set(it.product((-1, 0, 1), repeat=3))          # create seed
    self.grid = Grid(twins)
    self.create_entire_surface()


  def __repr__(self):
    return "fLakE :: twiNpLaNes:{}\tIterations:[{}]\tSeed:[3x3x3]".format(
      self.twins, self.iter)


# # ###########
#     BASIC   #
# #############
  def set_surface(self, site):
    """ Adds an empty site to the corresponding surface slot.
    """
    occupied_neighbours = self.real_neighbours(site, void=False)
    slot = len(occupied_neighbours)
    self.surface[slot].add(site)


  def create_entire_surface(self):
    """ Create a list of 12 sets (NN possibilities) and populates them.

    Iterates through every atom and populates each adjacent empty site to the
    surface list.
    """
    self.surface = [set() for _ in self.span]
    for atom in self.atoms:
      adjacent_voids = self.real_neighbours(atom, void=True)
      for nb in adjacent_voids:
        self.set_surface(nb)


  def geometry(self):
    """ Calculate geometry information about the Flake.

    Area is approximated as circle. Its radius is the distance from center to
    the furthest atom.
    """
    COORD = self.grid.coord
    POOL = self.atoms

    mxz = max(POOL, key=lambda i: i[2])[2]
    mnz = min(POOL, key=lambda i: i[2])[2]
    mxr = max(POOL, key=lambda i: i[0]**2 + i[1]**2)

    items = {
      'height': COORD((0, 0, mxz)).dist(COORD((0, 0, mnz))),
      'radius': COORD(mxr).dist(COORD((0, 0, mxr[2])))}
    items.update({
      'area': pi*items['radius']**2,
      'aspect_ratio': items['radius']/items['height']})

    for (k, i) in items.iteritems():
      setattr(self, k, i)
    return items


# # #################
#   NEIGHBORHOOD  #
# ##################
  def abs_neighbours(self, idx):
    """ Return a list of next neighbours of `i, j, k`
    """
    absNB = set(it.product((-1, 0, 1), repeat=3))     # create NN indices
    absNB.remove((0, 0, 0))
    pairs = (zip(idx, nn) for nn in absNB)
    absolutes = (tuple(sum(y) for y in x) for x in pairs)
    return absolutes


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
    indexed = list(self.abs_neighbours(atom))
    coordinates = (self.grid.coord(ai_site) for ai_site in indexed)
    diffs = (choice_vec.dist(ac_site) for ac_site in coordinates)
    associated = list(zip(indexed, diffs))
    nearest = [nb for nb, diff in associated if diff < DIFF_CAP]
    if void:
      return [nb for nb in nearest if nb not in self.atoms]
    else:
      return [nb for nb in nearest if nb in self.atoms]


############
#  GROWTH  #
############
  def grow(self, rounds=1, noise=10e-4):
    """ Transform a surface site into an atom.

    The site is chosen randomly from the highest populated surface slot. With a
    probability of `noise` the atoms will be chosen randomly from all available
    surface sites.
    """
    for r in range(rounds):
      if random() < noise:
        chosen = choice(list(it.chain.from_iterable(self.surface)))
        slot = next(i for i, x in enumerate(self.surface) if chosen in x)
        qprint("Random-Add: {at} in Slot: [{sl}]\nWe are @{it}".format(
          at=chosen, sl=slot, it=self.iter), quiet=Q)
      else:
        for slot in range(11, 0, -1):
          if self.surface[slot]:
            chosen = choice(tuple(self.surface[slot]))
            break
      self.put_atom(chosen, slot)
    self.geometry()
    return chosen

  def put_atom(self, at, slot):
    """ * remove atom from its slot
        * append it to atoms list
        * check which nb are free, with those:
          * check every slot if it contains the nb
          * remove it from there
          * add it to next higher slot (cause now he's got one nb more)
    """
    self.surface[slot].remove(at)
    self.atoms.add(at)
    empty_neighbours = self.real_neighbours(at, void=True)
    for each in empty_neighbours:
      for e_slot, lst in enumerate(self.surface):
        if each in lst:
          self.surface[e_slot].remove(each)
          try:
            self.surface[e_slot + 1].add(each)
          except IndexError:
            qprint("Filled a bubble...oO", quiet=True)
          break
      #TODO: proper else?
      else:
        self.surface[1].add(each)    # create new surface entry for new ones
    self.iter += 1


###########
#   PLOT  #
###########
  def plot(self, save=False, tag=''):
    """ The `color()` function appends values for colors to the `color_list`.
    This list must have the same length as `whole` to plot correctly.

    """
    surface_chain = list(it.chain.from_iterable(self.surface))
    whole = list(self.atoms) + surface_chain
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
      m.options.offscreen = True    # this should suppress output on screen
      if tag:                       # currently not working (bug in mayavi?)
        tag = str(tag) + '_'
      save_dir = self.daily_output()
      _time = self.date[1].rsplit('.')[0].replace(':', '-')
      fname = path.join(save_dir, 'Flake@' + _time + '_T' + str(self.twins) +
                        tag + '.png')
      m.savefig(fname, size=(1024, 768))
      m.close()
    else:
      m.show()

  def daily_output(self):
    self.date = arrow.now().isoformat().rsplit('T')
    today = self.date[0]
    output_dir = path.join('../output/', today)
    if not path.exists(output_dir):
      mkdir(output_dir)
    return output_dir


  def export(self, tag='flake'):
    """ Simplified export function adapted from 'io_mesh_xyz'.
    """
    raw_atoms = (('Au', tuple(self.grid.coord(at)))
                 for at in self.atoms if at)
    list_atoms = []
    counter = 0

    for each in raw_atoms:
      list_atoms.append(AtomsExport(*each))
      counter += 1

    save_dir = self.daily_output()
    fname = "{tag}._TP-{tp}_it-{it}_{tag}.xyz".format(
      tp=self.twins, it=len(self.atoms), tag=tag)
    filepath_xyz = path.join(save_dir, fname)
    with open(filepath_xyz, "w") as xyz_file_p:
      xyz_file_p.write("%d\n" % counter)
      xyz_file_p.write(str(self.geometry()))
      xyz_file_p.write("This is a XYZ file for Atomic Blender\n")

      for i, atom in enumerate(list_atoms):
          string = "%3s%15.5f%15.5f%15.5f\n" % (
                                        atom.element,
                                        atom.location[0],
                                        atom.location[1],
                                        atom.location[2])
          xyz_file_p.write(string)
