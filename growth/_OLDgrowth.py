#!/usr/bin/env python3
# encoding: utf-8              A flake growth simulation

from math import sqrt
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# ---------------
# -  The Flake  -
# ---------------
class Flake():
  """ The actual Flake we want to let grow. """
  def __init__(self):
    # lattice base vector, transposed to match [a, b, c] structure
    self.fcc_base = np.array(((1, 0, 0),
                           (1/sqrt(2), 1/sqrt(2), 0),
                           (1/2, 1/(2*sqrt(3)), sqrt(2/3))
                           )).T
    self.hcp_base = np.array(((1, 0, 0),
                           (1/sqrt(2), 1/sqrt(2), 0),
                           (1/2, 1/(2*sqrt(3)), sqrt(2/3))
                           )).T

  def plot(self, points=None):
    """ Plot method of the flake. Should point to lattice as default. """
    if points is None:
      points = self.fcc_base
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(*points, s=1500)
    plt.show(block=True)

  def permutator(self, seed, inplane=True):
    if inplane:
      dimsize = 2
    else:
      dimsize = 3
    types = it.combinations_with_replacement(seed, dimsize)
    perms = []
    for i in types:
      t = set(it.permutations(i))
      while t:
        perms.append(t.pop())
    return tuple(perms)

  def constructor(self, step_vector, basis='inplane'):
    """ Create x, y, z coordinates from the basis vector and translation
    number.
    """
    if basis == 'inplane':
      basis = self.fcc_base[:, :-1]
    elif basis == 'full':
        basis = self.fcc_base
    try:
      coords = np.dot(basis, step_vector)
    except:
      print("\nBase:\n", basis)
      print("\nStep:", step_vector)
      raise
    return coords.T

  def plane(self, z=0, size=2):
    """ Creates a plane of gold atoms.  `sites` create a list of all possible
    2-tuple combinations of integers in the range of `size`.
    """
    sites = self.permutator(range(1 - size, size))
    print("\nSITES: \n", sites)
    ret = np.zeros(2)
    for site in sites:
      co = self.constructor(site)
      ret = np.vstack((ret, co[:-1]))
    print("RET, COORDS:\n", ret, co)
    return ret[1:]

  def builder(self, indices, struct='fcc'):
    """ Build crystal as i*a + j*b + k*c with lattice vectors. """
    if struct == 'fcc':
      switch = 2
    elif struct == 'hcp':
      switch = 3
    else:
      raise NotImplementedError(
        struct + " is not valid.\nPlease choose 'fcc' or 'hcp'.")

    i, j, k = indices
    prototype = np.array((2*i + (j+k)%2,
                          sqrt(3)*(j + k%switch*1/3),
                          k*2*sqrt(6)/3
                          ))
    return prototype


  # def stack(height=5, size=3):
    # for plane in range(height):


# ----------
# -  Main  -
# ----------
def main():
  f = Flake()
  pp = f.plane(size=5)
  print("planes;", pp)
  f.plot(pp)


if __name__ == '__main__':
  main()
  if False:
    Axes3D
