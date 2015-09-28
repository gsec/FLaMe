#!/usr/bin/env python3
# encoding: utf-8              A flake growth simulation

from math import sqrt
import numpy as np
from numpy import matrix
import itertools as it
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Flake():
  """ The actual Flake we want to let grow"""
  def __init__(self):
    # self.lattice  = np.zeros(5 * self.size**3).reshape(
                    # self.size,      # x'
                    # self.size,      # y'
                    # self.size,      # z'
                    # 5)              # (x, y, z, t, u)

    self.basis_vector = matrix(((1, 0, 0),
                                (1/sqrt(2), 1/sqrt(2), 0),
                                (1/2, 1/(2*sqrt(3)), sqrt(2/3))))

  def coordinates(self, ordinals, vector_set=None):
    """ Create x, y, z coordinates from the basis vector and translation
    number.
    """
    if vector_set is None:
      vector_set = self.basis_vector
    if len(ordinals) == 2:
      vector_set = vector_set[:2]
    coords = matrix(ordinals)*matrix(vector_set)
    return np.array(coords)

  def plot(self, points=None):
    if points is None:
      points = self.basis_vector
    points = np.asarray(points).T

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    print("\npoints", points)
    ax.scatter(*points, s=1500)
    plt.show()

  def plane(self, z=0, size=2):
    """ Creates a plane of gold atoms.
    """
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


if False:
  Axes3D

if __name__ == '__main__':
  main()
