#!/usr/bin/env python3
# encoding: utf-8
"""                              FLAKE GROWTH SIMULATION
                                 ~~~~~~~~~~~~~~~~~~~~~~~

        :copyright: (c) 2015     Guilherme Stein
                                 University of Würzburg
                                 <guilherme.stein@physik.uni-wuerzburg.de>
"""

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
  def __init__(self, size: int=7, twins=tuple()):
    self.Vector = namedtuple('Vector', ['x', 'y', 'z'])
    self.grid_list = [[[[0, 0] for _ in range(size)]
                       for _ in range(size)]
                      for _ in range(size)]

    self.size = size
    self.layer_permutations = self.layer_gen(*twins)

  def grid(self, i, j, k, val=None):
    """ Interface to acess the `grid_list` containing information about each atom.
    """
    point = self.grid_list[i][j][k]

    if val is None:
      return point[0]
    elif val == 'energy':
      return point[1]
    elif val == 'full':
      return point
    elif val in (True, False):
      self.grid_list[i][j][k] = val
      return point
    elif isinstance(val, (list, tuple)) and len(val) == 2:
      self.grid_list[i][j][k] = val
    else:
      raise TypeError("Unrecognized type of `val`.\nMust be either `keyword`, "
                      "`boolean` or an iterable with length two.")

  def coord(self, i, j, k) -> 'Vector':
    """ Build crystal as i*a + j*b + k*c with lattice vectors.

    `perms` is the permutation (0, 1 or 2) of the layer displacement according to the
    fcc-stacking and the twin plane configuration (every twin plane inverts the
    permutation order).  """
    perms = self.layer_permutations[k]
    prototype = self.Vector(2*i + (j+k) % 2,
                            sqrt(3)*(j + perms * 1/3),
                            k*2*sqrt(6)/3)
    return prototype

  def layer_gen(self, *twins) -> list:
    """ Create a z-list representing the permutation of the layer.
    """
    L = []
    sign = 1
    counter = 0

    for layer in range(self.size):
      L.append(counter % 3)
      if layer in twins:
        sign = -1*sign
      counter += sign
    return L

  def plot(self, points: list, color: list=None):
    """ Plot method of the flake. """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    pts = list(zip(*points))
    if not color:
      # color = pts[2]
      color = []
    ax.scatter(*pts, s=1500, c=color)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

  def nn_gen(self, i, j, k, stack_list=None):
    """ Return the combination of next neighbours diff vectors, depending on the upper and
    lower plane indices.

    TP: twin plane
    U(D): Up/Down
    N: normal plane
    """
    IN_PLANE = [(1, 0, 0), (-1, 0, 0), (0, 1, 0),
                (0, -1, 0), (-1, -1, 0), (-1, 1, 0)]
    N_U_N = [(0, 0, 1), (-1, 0, 1), (0, 1, 1)]
    N_D_N = [(0, 0, -1), (-1, 0, -1), (0, 1, -1)]
    TP_D_N = [(0, 0, -1), (1, 0, -1), (0, -1, -1)]    # same as TP_D_TP
    TP_U_N = [(0, 0, 1), (1, 0, 1), (0, 1, 1)]        # same as TP_U_TP
    N_U_TP = [(0, 0, 1), (1, 0, 1), (0, -1, 1)]
    N_D_TP = [(0, 0, -1), (-1, 0, -1), (0, 1, -1)]

    if not stack_list:
      stack_list = self.layer_permutations
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

  def set_neighbours(self, i, j, k, val=1):
    for n in self.neighbours(i, j, k):
      self.grid(*n, val=val)

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

  def all_points(self):
    coords = []
    colors = []
    for idx in self.permutator(range(self.size)):
      _c = self.coord(*idx)
      coords.append(_c)
      colors.append(self.grid(*idx, val='energy') + 0.1 * idx[2])
    return coords, colors


# ----------
# -  main  -
# ----------
def main():
  f = Flake(size=3)
  coords = []
  cols = []
  if Axes3D:
    special_one = (1, 1, 1)
    f.grid(*special_one, val=(1, 5))
    f.set_neighbours(*special_one, val=(1, 3))
    print('\nRendering:')
    for idx in f.permutator(range(f.size)):
      # if f.grid(*idx):
        _c = f.coord(*idx)
        coords.append(_c)
        cols.append(f.grid(*idx, val='energy') + 0.1 * idx[2])
        # print(idx, _c, sep='\t')
    f.plot(coords, cols)


if __name__ == '__main__':
  main()
