#!/usr/bin/env python2
# encoding: utf-8
"""                           FLagS - a FLakeLatticeGrowthSimulator
                              <guilherme.stein@physik.uni-wuerzburg.de>

"""

# from flame.growth import Flake
# from random import random
# from os import getcwd, path
# import tables as tb
# import logging
# import yaml
# import sys

# LOG CONFIG
# logging.basicConfig(level=logging.INFO)
# logger = logging.getLogger(__name__)

from __future__ import print_function, division, generators
from flame.simulation import Stage
from os import environ, path

try:
    SIM_DIR_GLOBAL = path.join(environ['THESIS_PATH'], 'output', 'sim')
except KeyError as e:
    logger.warn("No {} specified. Please specify path manually.".format(e))
    SIM_DIR_GLOBAL = ''
    while not path.isdir(SIM_DIR_GLOBAL):
        SIM_DIR_GLOBAL = input("Top level simulation directory: ")
