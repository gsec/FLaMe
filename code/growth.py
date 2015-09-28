#!/usr/bin/env python3
# encoding: utf-8
# A flake growth simulation
# from collections import namedtuple
from math import sqrt
import numpy as np


class Flake():
  """ The actual Flake we want to let grow"""
  def __init__(self):
    # self.Vector = namedtuple('Vectors', ['a', 'b', 'c'])
    self.Vector = np.asarray
    self.basis_vector = self.Vector((
                                   (1, 0, 0),
                                   (1/sqrt(2), 1/sqrt(2), 0),
                                   (1/2, 1/(2*sqrt(3)), sqrt(2/3))
                                   ))

  def coordinates(self, ordinals):
    """ Create x, y, z coordinates from the basis vector and translation
    number.
    """
    res = []
    for idx, item in enumerate(self.Vector(ordinals)):
      ret = item * self.basis_vector[idx]
      res.append(ret)

    # res = self.basis_vector * self.Vector(ordinals)
    print("Multi", res)
    # return (res).sum(0)
