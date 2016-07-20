#!/usr/bin/env python2
#coding: utf-8
#                     The test module for the flake analysis
from __future__ import print_function, division, generators
import unittest
import analysis as A
from os.path import isdir




class TestReadout(unittest.TestCase):
    def SetUp(self):
        pass

    def test_output_dir_existence(self):
        self.assertTrue(isdir(A.SIM_DIR_GLOBAL))

if __name__ == '__main__':
    unittest.main()
