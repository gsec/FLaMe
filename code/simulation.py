#!/usr/bin/env python2
# encoding: utf-8

# Define the general Structure of the Simulations (here we use this as a
# synonym to experiments). We will create a hierarchy in the HDF data format and
# group those elements structure for the experiments we want to run with our basic flake
# growth module

from __future__ import print_function, division, generators
from growth import Flake
from random import random
from os import path, environ
import logging
import arrow
import h5py

# LOG CONFIG
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

CFG = {
    'SAMPLE_SIZE': 5,
    'SNAPSHOT_INTERVAL': 1000,
    'TOTAL_SIZE': 10000,
    'UNIT_NAME': 'Ycaruz',
    'UNIT_SPAN': 3,
    'OUTPUT_DIR': path.join(environ['THESIS_PATH'], 'output/sim')
    }
# Configuration variables:
SAMPLE_SIZE = 5             # number of flakes grown the same to average
SNAPSHOT_INTERVAL = 1000    # number of atoms to grow between flake snapshots
TOTAL_SIZE = 10000          # total number of atoms at the end of growth
UNIT_NAME = 'Ycaruz'        # name of the simulation
UNIT_SPAN = 3               # span over which chosen parameter is varied
OUTPUT_DIR = path.join(environ['THESIS_PATH'], 'output/sim')


# Create a new Dataset
def create(name):
    """ Spawn HDF file and create basic structure for the experiment.
    """
    identifier = str(name) + '_' + hex(hash(random()))
    logger.info('\nSIMULATION STARTED >>> {}'.format(identifier))
    for k, v in CFG.items():
        logger.info('{}: {}'.format(k, v))
    fname = path.join(OUTPUT_DIR, identifier + '.hdf5')
    f = h5py.File(fname, 'x')

    mapping = lambda x: (-x, x + 1)

    for num in range(UNIT_SPAN):
        tp_set = mapping(num)

        param_grp = f.create_group('param=' + str(num))
        for smp in range(SAMPLE_SIZE):
            time = '@' + arrow.now().format('YYYY\'MM\'DD HH:mm:ss')
            logger.info(time + ' Sampling... {}/{}'.format(smp+1, SAMPLE_SIZE))
            header, gdata = run(tp_set)
            sample_grp = param_grp.create_group('sample{:04}'.format(smp))
            sample_grp.create_dataset('flake', data=gdata)
            if num == 0 and smp == 0:
                f.attrs['header'] = str(header)
    f.close()
    logger.info('SIMULATION ENDED >>> {}'.format(identifier))


def run(tp, total_size=TOTAL_SIZE, samples=SAMPLE_SIZE, intervall=1000):
    myFlake = Flake(*tp)
    div, remainder = total_size // intervall, total_size % intervall

    header = myFlake.geometry().keys()
    data = []

    logger.info('Twinplanes: {} \tTotal Size:{} \tSeed: {} \tTemp:[{}]'.format(
      myFlake.twins, total_size, myFlake.seed_shape, myFlake.temp))

    for it in range(div):           # for all whole steps
        myFlake.grow(intervall)
        data.append(myFlake.geometry().values())
    if remainder:                   # for the remaining ones
        myFlake.grow(remainder)
        data.append(myFlake.geometry().values())

    return (header, data)
