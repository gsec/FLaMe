#!/usr/bin/env python2
## encoding: utf-8
"""                       FLAKE GROWTH SIMULATION
                          ~~~~~~~~~~~~~~~~~~~~~~~

                          Guilherme Stein © 2015
                          University of Würzburg
                          <guilherme.stein@physik.uni-wuerzburg.de>
"""
from __future__ import print_function, division, generators
import arrow
import itertools as it
from os import path, mkdir
from collections import deque, OrderedDict
from math import pi
from helper import *
from mayavi import mlab as m
from numpy.random import choice as weighted_choice
from random import choice, random

Q = False          # Set verbosity, Q is quiet


class Flake(object):
  """ FLAKE creation class.
  """
  def __init__(self, *twins, **kwargs):
    self.twins = twins
    self.span = range(12)                 # take all 12 NN as possibilities

    self.seed_shape = kwargs.get('seed', 'sphere')
    self.trail_length = kwargs.get('trail', 10)
    self.temp = kwargs.get('temp', 10)
    self.iter, self.atoms = self.seed(self.seed_shape)
    self.trail = deque(maxlen=self.trail_length)   # specify length of trail

    self.grid = Grid(twins)
    self.create_entire_surface()
    self.colors = self.color_init()


  def __repr__(self):
    return ("fLakE :: twiNpLaNes:{}\tIterations:[{}]\t Seed:[{}]\tTemperature:"
            "[{}K]".format(self.twins, self.iter, self.seed_shape, self.temp))


# # ###########
#     BASIC   #
# #############
  def seed(self, shape='point'):
    if shape == 'point':
      sites = set(((0, 0, 0),))
    elif shape == 'cube':
      sites = set(it.product((-1, 0, 1), repeat=3))
    elif shape == 'bigcube':
      sites = set(it.product((-2, -1, 0, 1, 2), repeat=3))
    elif shape == 'sphere':
      sites = set(((0, -1, 1), (0, 0, 1), (-1, 0, 1), (-1, -1, -1), (0, -1, -1),
                  (0, 0, -1), (-1, -1, 0), (-1, 0, 0), (-1, 1, 0), (0, -1, 0),
                  (0, 1, 0), (1, 0, 0), (0, 0, 0)))
    elif shape == 'plane':
      sites = set(((1, 0, 0), (1, 2, 0), (-2, 1, 0), (-2, 0, 0), (1, -1, 0),
                  (0, 1, 0), (-2, 2, 0), (-1, 0, 0), (-2, -2, 0), (0, -1, 0),
                  (1, 1, 0), (1, -2, 0), (0, -2, 0), (0, 2, 0), (2, 0, 0),
                  (-1, -2, 0), (-1, 1, 0), (-1, 2, 0), (2, -1, 0), (-1, -1, 0),
                  (2, 2, 0), (0, 0, 0), (2, 1, 0), (-2, -1, 0), (2, -2, 0)))
    return len(sites), sites

  def sites(self):
    return OrderedDict(('Slot [{}]'.format(i), len(x)) for (i, x) in
                     enumerate(self.surface))

  def color_init(self):
    d = dict((k, 15) for k in self.atoms)
    d.update([(y, k) for k, x in enumerate(self.surface) for y in x])
    return d


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

    items = {                                         # +1 correct for border
      'height': COORD((0, 0, mxz + 1)).dist(COORD((0, 0, mnz))),
      'radius': COORD(mxr).dist(COORD((0, 0, mxr[2])))}
    items.update({
      'area': pi*items['radius']**2,
      'aspect_ratio': 2*items['radius']/items['height']})

    for (k, i) in items.iteritems():
      setattr(self, k, i)
    return items


# # #################
#   NEIGHBORHOOD  #
# ##################
  def abs_neighbours(self, idx):
    """ Return a list of next neighbours of `i, j, k`
    """
    absNB = set(it.product((-2, -1, 0, 1, 2), repeat=3))     # create NN indices
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
    DIFF_CAP = 2.1
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


#################
#  PROBABILITY  #
#################
  def prob(self):
    weights = []
    factor = 5.
    func = lambda x: x**((300 - self.temp)/factor) * len(self.surface[x])
    norm = sum(func(idx) for idx, item in enumerate(self.surface))
    for slot in self.span:
      p = func(slot) / norm
      weights.append(p)
    assert(1 - sum(weights) < 1e-12)  # Probability check
    return weights


  def cave(self):
    temp_atoms = self.atoms.copy()
    ret = []
    for atom in temp_atoms:
      if not self.real_neighbours(atom, void=True):
        self.atoms.remove(atom)
        ret.append(atom)
    return ret


############
#  GROWTH  #
############
  def grow(self, rounds=1):
    for r in range(rounds):
      sp = self.prob()
      slot = weighted_choice(self.span, p=sp)
      chosen = choice(tuple(self.surface[slot]))
      self.put_atom(chosen, slot)


  def det_grow(self, rounds=1, noise=0, cap=1):
    """ Transform a surface site into an atom.

    The site is chosen randomly from the highest populated surface slot. With a
    probability of `noise` the atoms will be chosen randomly from all available
    surface sites.
    """
    for r in range(rounds):
      if noise:
        if random() < noise:
          chosen, slot = self.rand_grow('NoiZ')
      else:
        for slot in range(11, cap-1, -1):
          if self.surface[slot]:
            chosen = choice(tuple(self.surface[slot]))
            break
        else:
          print("Flake has no sites with {} or more free bindings. STOPP at {}"
                "th growth step.\n".format(cap, r))
          ans = raw_input("Continue with single growth? [y/n]\t")
          if ans in 'Yy':
            self.det_grow()
            continue
          break
      self.put_atom(chosen, slot)
    return self.trail

  def rand_grow(self, text):
    chosen = choice(list(it.chain.from_iterable(self.surface)))
    slot = next(i for i, x in enumerate(self.surface) if chosen in x)
    msg = "[{}] add {} in Slot: [{}]\nWe are @ {} iterations.\n".format(text,
      chosen, slot, self.iter)
    qprint(msg, quiet=Q)
    return chosen, slot

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
    self.colors.update(((at, 11),))
    if len(self.trail) == self.trail_length:
      old = self.trail.pop()
      self.colors.update(((old, 15),))
    self.trail.appendleft(at)        # create list of latest additions

    empty_neighbours = self.real_neighbours(at, void=True)
    for each in empty_neighbours:
      for e_slot, lst in enumerate(self.surface):
        if each in lst:
          self.surface[e_slot].remove(each)
          try:
            self.surface[e_slot + 1].add(each)
            self.colors.update(((each, e_slot + 1),))
          except IndexError:
            # self.colors.remove(((each, e_slot + 1),))
            self.colors.pop(each)
            qprint("Filled a bubble...oO", quiet=Q)
          break
      else:
        self.surface[1].add(each)    # create new surface entry for new ones
    self.iter += 1


#############
#   OUTPUT  #
#############
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


  def plot(self, save=False, tag='', pipeline=False):
    """ The `color()` function appends values for colors to the `color_list`.
    This list must have the same length as `whole` to plot correctly.

    """
    if pipeline:
      co = self.grid.coord
      clist = zip(*[(co(idx).x, co(idx).y, co(idx).z, c) for (idx, c) in self.colors.items()])
      # clist = [(co(idx).x, co(idx).y, co(idx).z, c) for (idx, c) in self.colors.items()]
      return clist

    surface_chain = list(it.chain.from_iterable(self.surface))
    whole = list(self.atoms) + surface_chain
    x, y, z = list(zip(*(self.grid.coord(site) for site in whole)))

    color_list = []
    for each in self.atoms:
      if each in self.trail:
        val = 10
      else:
        val = 15
      color_list.append(val)

    for (idx, surf) in enumerate(self.surface):
      for each in surf:
        color_list.append(1.*idx)

    # m.clf()
    m.points3d(x, y, z, color_list, colormap="gist_ncar",
                         scale_factor=0.1, vmin=0, vmax=15)
    if save == 1:
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
      return x, y, z, color_list
      # m.show()
