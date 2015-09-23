#!/usr/bin/env python3
# encoding: utf-8
# A flake growth simulation
from collections import namedtuple


class Flake():
  """ The actual Flake we want to let grow"""
  def __init__(self):
    vectors = namedtuple('Vectors', ['a', 'b', 'c'])
    self.vectors = vectors((0, 0, 0), (0, 0, 0), (0, 0, 0))
