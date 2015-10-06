#!/usr/bin/env python3
# encoding: utf-8              A flake growth simulation
# pylint: disable=W0611

from math import sqrt
# import numpy as np
import itertools as it
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import namedtuple


# ---------------
# -  The Flake  -
# ---------------
class Flake():
  """ The actual Flake we want to let grow. """
  def __init__(self, size=5):
    self.size = size
    self.Vector = namedtuple('Vector', ['x', 'y', 'z'])
    self.LatticeBase = namedtuple('Base', ['a', 'b', 'c'])
    a = self.Vector(1, 0, 0)
    b = self.Vector(1/sqrt(2), 1/sqrt(2), 0)
    c = self.Vector(1/2, 1/(2*sqrt(3)), sqrt(2/3))
    self.fcc_base = self.LatticeBase(a, b, c)
    y_hcp = self.coord(1, 2, 3, 'hcp').y
    y_fcc = self.coord(1, 2, 3, 'fcc').y
    self.lattice_delta = self.Vector(0, y_hcp - y_fcc, 0)
    self.grid_list = [[[False for _ in range(size)]
                       for _ in range(size)]
                      for _ in range(size)]
    self.layers = [0 for _ in range(size)]

  def grid(self, i, j, k, val=None):
    if val is None:
      return self.grid_list[i][j][k]
    elif val in (True, False):
      self.grid_list[i][j][k] = bool(val)
      return self.grid_list[i][j][k]
    else:
      raise TypeError("Value must be boolean!")

  def coord(self, i, j, k, struct=None):
    """ Build crystal as i*a + j*b + k*c with lattice vectors. """
    if not struct:
      struct = self.layers[j]
    elif struct in ('fcc', 0):
      struct = 0
    elif struct in ('hcp', 1):
      struct = 1
    else:
      raise KeyError("The structure `" + str(struct)  + "` is not defined")
    # This is the modulo base depending on the lattice structure:
    switch = {0: 2, 1: 3}[struct]
    prototype = self.Vector(2*i + (j+k) % 2,
                            sqrt(3)*(j + k % switch * 1/3),
                            k*2*sqrt(6)/3)
    return prototype

  def plot(self, points):
    """ Plot method of the flake. """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(*points, s=1500)
    plt.show()

  def permutator(self, seed):
    """ Creates all possible permutations of length three of all given objects
    in `seed`. """
    types = it.combinations_with_replacement(seed, 3)
    perms = []
    for i in types:
      t = set(it.permutations(i))
      while t:
        perms.append(t.pop())
    return perms


def main():
  SIZE = 7
  f = Flake(SIZE)
  coords = []
  for idx in f.permutator(range(-SIZE//2, SIZE//2)):
    print(idx)
    _c = f.coord(*idx, struct=1)
    coords.append(_c)
    print(_c)
  print(coords)
  f.plot(zip(*coords))


if __name__ == '__main__':
  main()
if Axes3D:
  pass
