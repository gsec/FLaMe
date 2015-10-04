#!/usr/bin/env python3
# encoding: utf-8              A flake growth simulation
# pylint: disable=W0611

from math import sqrt
# import numpy as np
import itertools as it
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from collections import namedtuple


# ---------------
# -  The Flake  -
# ---------------
class Flake():
  """ The actual Flake we want to let grow. """
  def __init__(self, size=5):
    self.Vector = namedtuple('Vector', ['x', 'y', 'z'])
    self.LatticeBase = namedtuple('Base', ['a', 'b', 'c'])

    a = self.Vector(1, 0, 0)
    b = self.Vector(1/sqrt(2), 1/sqrt(2), 0)
    c = self.Vector(1/2, 1/(2*sqrt(3)), sqrt(2/3))
    self.fcc_base = self.LatticeBase(a, b, c)

    self.grid_list = [[[False for _ in range(size)]
                       for _ in range(size)]
                      for _ in range(size)]

  def grid(self, i, j, k, val=None):
    if val is None:
      return self.grid_list[i][j][k]
    if val not in (True, False):
      raise TypeError("Value must be boolean!")
    self.grid_list[i][j][k] = bool(val)

  def coord(self, i, j, k, struct='fcc'):
    """ Build crystal as i*a + j*b + k*c with lattice vectors. """
    switch = {'fcc': 2, 'hcp': 3}[struct]
    prototype = self.Vector(2*i + (j+k) % 2,
                            sqrt(3)*(j + k % switch * 1/3),
                            k*2*sqrt(6)/3)
    return prototype

  def permutator(self, seed):
    types = it.combinations_with_replacement(seed, 3)
    perms = []
    for i in types:
      t = set(it.permutations(i))
      while t:
        perms.append(t.pop())
    return perms
