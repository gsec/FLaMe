#!/usr/bin/env python3
# encoding: utf-8              A flake growth simulation

from math import sqrt
import itertools as it
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import namedtuple


# ---------------
# -  The Flake  -
# ---------------
class Flake():
  """ The actual Flake we want to let grow. """
  def __init__(self, size=7, twins=None):
    self.size = size
    self.Vector = namedtuple('Vector', ['x', 'y', 'z'])
    self.LatticeBase = namedtuple('Base', ['a', 'b', 'c'])
    a = self.Vector(1, 0, 0)
    b = self.Vector(1/sqrt(2), 1/sqrt(2), 0)
    c = self.Vector(1/2, 1/(2*sqrt(3)), sqrt(2/3))
    self.fcc_base = self.LatticeBase(a, b, c)
    self.grid_list = [[[False for _ in range(size)]
                       for _ in range(size)]
                      for _ in range(size)]
    self.layers = self.layer_generator(twins)

  def grid(self, i, j, k, val=None):
    if val is None:
      return self.grid_list[i][j][k]
    elif val in (True, False):
      self.grid_list[i][j][k] = bool(val)
      return self.grid_list[i][j][k]
    else:
      raise TypeError("Value must be boolean!")

  def coord(self, i, j, k, twin=None):
    """ Build crystal as i*a + j*b + k*c with lattice vectors.
    `twin` is the variable that is increased by one for every twin plane we
    introduce. All layers _after_ the twin (greater positive and negative
    values) are mirrored. """
    if twin is None:
      twin = self.layers[k]
    layer = ((-1)**twin * k - twin) % 3
    # -1**t is to invert the order, -t to displace the index properly
    prototype = self.Vector(2*i + (j+k) % 2,
                            sqrt(3)*(j + layer * 1/3),
                            k*2*sqrt(6)/3)
    return prototype

  def layer_generator(self, twins=None):
    layers = [0 for _ in range(self.size)]
    if twins:
      if isinstance(twins, int):
        twins = [twins]
      if isinstance(twins, tuple):
        twins = list(twins)
      for i, v in enumerate(layers):
        layers[i] = len([t for t in twins if i > t])
    return layers

  def plot(self, points):
    """ Plot method of the flake. """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    points = list(points)
    ax.scatter(*points, s=1500, c=points[2])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.autoscale_view(None, False, False)
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
  SIZE = 3
  f = Flake(SIZE, twins=1)
  coords = []
  for idx in f.permutator(range(SIZE)):
    print(idx)
    _c = f.coord(*idx)
    coords.append(_c)
    print(_c)
  print(coords)
  f.plot(zip(*coords))


if __name__ == '__main__':
  main()
