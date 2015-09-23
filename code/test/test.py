#!/usr/bin/env python3
#coding: utf-8
# The test module for our program

import unittest
import growth


class TestBasicLattice(unittest.TestCase):
  """ the first test."""

  def setUp(self):
    """ Instantiate the class."""
    self.Flake = growth.Flake()

  def test_vectors(self):
    self.assertTrue(len(self.Flake.vectors) == 3)
    for vector in self.Flake.vectors:
      self.assertTrue(vector)

if __name__ == '__main__':
  unittest.main()
