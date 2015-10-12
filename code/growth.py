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
  def __init__(self, size: int=7, twins=None):
    self.size = size
    self.Vector = namedtuple('Vector', ['x', 'y', 'z'])
    # self.LatticeBase = namedtuple('Base', ['a', 'b', 'c'])
    # a = self.Vector(1, 0, 0)
    # b = self.Vector(1/sqrt(2), 1/sqrt(2), 0)
    # c = self.Vector(1/2, 1/(2*sqrt(3)), sqrt(2/3))
    # self.fcc_base = self.LatticeBase(a, b, c)
    self.grid_list = [[[False for _ in range(size)]
                       for _ in range(size)]
                      for _ in range(size)]
    self.layers = self.layer_generator(twins)

  def grid(self, i, j, k, val=None) -> bool:
    if val is None:
      return self.grid_list[i][j][k]
    elif val in (True, False):
      self.grid_list[i][j][k] = bool(val)
      return self.grid_list[i][j][k]
    else:
      raise TypeError("Value must be boolean!")

  def coord(self, i, j, k, twin=None) -> 'Vector':
    """ Build crystal as i*a + j*b + k*c with lattice vectors.
    `twin` is the variable that is increased by one for every twin plane
    we introduce. All layers _after_ the twin (greater positive and negative
    values) are mirrored. """
    if twin is None:
      twin = self.layers[k]
    layer = ((-1)**twin * k - twin) % 3
    # -1**t is to invert the order, -t to displace the index properly
    prototype = self.Vector(2*i + (j+k) % 2,
                            sqrt(3)*(j + layer * 1/3),
                            k*2*sqrt(6)/3)
    return prototype

  def layer_generator(self, twins=None) -> list:
    """ Create a z-list representing the number of twin-planes until up to this
    layer. """
    layers = [0 for _ in range(self.size)]
    if twins:
      if isinstance(twins, int):
        twins = [twins]
      if isinstance(twins, tuple):
        twins = list(twins)
      for i, v in enumerate(layers):
        layers[i] = len([t for t in twins if i > t])
    return layers

  def plot(self, points: 'Vector'):
    """ Plot method of the flake. """
    points = list(points)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    pts = zip(*points)

    ax.scatter(*pts, s=1500, c=points[2])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.autoscale_view(None, False, False)
    plt.show()

  def nn_gen(self, i, j, k, stack_list=None):
    """ Return the combination of next neighbours diff vectors, depending on
    the upper and lower plane indices.
    TP: twin plane; U(D): Up/Down; N: normal plane
    """
    IN_PLANE = [(1, 0, 0), (-1, 0, 0), (0, 1, 0),
                (0, -1, 0), (1, -1, 0), (-1, 1, 0)]
    N_U_N = [(0, 0, 1), (-1, 0, 1), (0, -1, 1)]
    N_D_N = [(0, 0, -1), (-1, 0, -1), (0, 1, -1)]
    TP_D_N = [(0, 0, -1), (1, 0, -1), (0, -1, -1)]    # same as TP_D_TP
    TP_U_N = [(0, 0, 1), (1, 0, 1), (0, 1, 1)]        # same as TP_U_TP
    N_U_TP = [(0, 0, 1), (1, 0, 1), (0, -1, 1)]
    N_D_TP = [(0, 0, -1), (1, 0, -1), (0, 1, -1)]

    if not stack_list:
      stack_list = self.layers
    try:
      tp_next = stack_list[k+1] - stack_list[k]
    except IndexError:
      print("Upper lattice border reached")
      tp_next = 0
    is_tp = stack_list[k] - stack_list[k-1]
    try:
      tp_prev = stack_list[k-1] - stack_list[k-2]
    except IndexError:
      print("Lower lattice border reached")
      tp_prev = 0

    neighbours = IN_PLANE
    if is_tp:
      neighbours.extend(TP_U_N)
      neighbours.extend(TP_D_N)
    if not is_tp:
      if tp_next:
        neighbours.extend(N_U_TP)
        if tp_prev:
          neighbours.extend(N_D_TP)
        elif not tp_prev:
          neighbours.extend(N_D_N)
      elif not tp_next:
        neighbours.extend(N_U_N)
        if tp_prev:
          neighbours.extend(N_D_TP)
        elif not tp_prev:
          neighbours.extend(N_D_N)
    return neighbours

  def neighbours(self, i, j, k):
    """ Return a list of next neighbours of `i, j, k` """
    site = (i, j, k)
    nn_vec = self.nn_gen(i, j, k)
    pairs = [zip(site, nn) for nn in nn_vec]
    return [tuple(sum(y) for y in x) for x in pairs]

  def set_neighbours(self, i, j, k):
    for n in self.neighbours(i, j, k):
      self.grid(*n, val=True)

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


# ------------------
# -  Global Stuff  -
# ------------------
# class LayerError(Exception):
  # pass


if __name__ == '__main__':
  main()
