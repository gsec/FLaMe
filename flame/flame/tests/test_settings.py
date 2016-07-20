#!/usr/bin/env python2
#coding: utf-8
#                     The test module for the settings
from __future__ import print_function, division, generators
import unittest
from flame import growth
from flame import settings as st
from flame.helper import Vector
from os import environ, path


class TestEnvironment(unittest.TestCase):
    """ Check environ variables and directories. """
    def SetUp(self):
        """ Clear environment for testing purposes.
        """
        self.original_environ = environ
        environ.clear()

        self.THESIS_PATH  = 'repos/global/thesis'
        self.FLAME_OUTPUT = 'repos/global/thesis/output'
        self.GROW_OUTPUT  = 'repos/global/thesis/output/grow'
        self.SIM_OUTPUT   = 'repos/global/thesis/output/sim'

    def test_output_dir(self):
        growth.FLAME_OUTPUT = ''
        with self.assertRaises(KeyError):
            print(growth.FLAME_OUTPUT)
        # self.assert(path.isdir())



