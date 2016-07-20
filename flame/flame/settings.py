""" FLaMe settings
"""
from os import environ, path
from logging import getLogger


logger = getLogger(__name__)
global FLAME_OUTPUT, GROW_OUTPUT, DIFF_CAP, PICKLE_EXT, HDF_EXT

try:
    FLAME_OUTPUT = environ['FLAME_OUTPUT']
except KeyError as e:
    logger.warn("No {} specified. Please specify path manually.".format(e))
    FLAME_OUTPUT = ''
    while not path.isdir(FLAME_OUTPUT):
        FLAME_OUTPUT = input("Top level output directory: ")
logger.info("\nFLaMe output path set to: {}".format(FLAME_OUTPUT))


GROW_OUTPUT = path.join(FLAME_OUTPUT, 'grow')
DIFF_CAP = 2.1
PICKLE_EXT = '.flm'
HDF_EXT = '.hdf'
