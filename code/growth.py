#!/usr/bin/env python3
# encoding: utf-8
# A flake growth simulation
from collections import namedtuple
from math import sqrt


class Flake():
  """ The actual Flake we want to let grow"""
  def __init__(self):
    self.Vector = namedtuple('Vectors', ['a', 'b', 'c'])
    self.basis_vector = self.Vector(
                                   (1, 0, 0),
                                   (1/sqrt(2), 1/sqrt(2), 0),
                                   (1/2, 1/(2*sqrt(3)), sqrt(2/3))
                                   )
