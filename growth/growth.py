#!/usr/bin/env python3
# encoding: utf-8              A flake growth simulation

from math import sqrt
import numpy as np
from numpy import array, matrix
import itertools as it
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Flake():
  """ The actual Flake we want to let grow. """
  def __init__(self):
    # lattice base vector, transposed to match [a, b, c] structure
    self.fcc_base = array(((1, 0, 0),
                           (1/sqrt(2), 1/sqrt(2), 0),
                           (1/2, 1/(2*sqrt(3)), sqrt(2/3))
                          )).T

  def coordinates(self, step_vector, base_vectors=None):
    """ Create x, y, z coordinates from the basis vector and translation
    number.
    """
    if base_vectors is None:
      base_vectors = self.fcc_base
    if len(step_vector) == 2:
      base_vectors = base_vectors[:,:-1]
      print("A Short one", base_vectors)
    # Transpose, so we have each component in proper slot of return vector
    # import pdb; pdb.set_trace()
    coords = matrix(base_vectors) * matrix(step_vector).T
    return np.array(coords)

  def plot(self, points=None):
    if points is None:
      points = self.fcc_base
    points = np.asarray(points).T

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    print("\npoints", points)
    ax.scatter(*points, s=1500)
    plt.show()

  def plane(self, z=0, size=2):
    """ Creates a plane of gold atoms. """
    sites = list(it.combinations_with_replacement(range(- size, size), 2))
    # a = np.array(())
    print("site combin", sites)
    a = np.zeros(2)
    # print(a.ndim)
    for site in sites:
      print("before", a)
      # a = np.concatenate((a[:,np.newaxis], site[:,np.newaxis]), axis=1)
      # a = np.concatenate((a, site), axis=0)
      co = self.coordinates(site)
      # import pdb; pdb.set_trace()
      # print("len", len(a), len(co))
      # print(co, type(co))
      # print(co[:2])
      a = np.vstack((a, co[:,:-1]))
      print("after", a)
    return a[1:]
    # return np.asarray(self.coordinates(site) for site in sites)


def main():
  f = Flake()
  # for atom in f.plane():
    # co = f.coordinates(atom)
  # f.plot(f.plane())
  # plt.show()
  pp = f.plane(size=5)
  print("planes;", pp)
  f.plot(pp)


if __name__ == '__main__':
  main()
  if False:
    Axes3D
