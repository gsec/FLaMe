#!/usr/bin/env python2
# encoding: utf-8
"""                           FLagS - a FLakeLatticeGrowthSimulator

                              <guilherme.stein@physik.uni-wuerzburg.de>
"""

# Define the general Structure of the Simulations (here we use this as a
# synonym to experiments). We will create a hierarchy in the HDF data format and
# group those elements structure for the experiments we want to run with our basic flake
# growth module

from __future__ import print_function, division, generators
from growth import Flake
from random import random
from os import getcwd, path
import argparse
import logging, arrow, h5py, yaml

# LOG CONFIG
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def read_params():
    """ Ensures the parameters are read in from a YAML file.

    A simulation is defined through its configuration file. Here we use YAML
    since it is human- and machine-readable.
    """
    def get_config():
        try:
            with open('tparams.yaml', 'r') as file_handler:
                params = yaml.load(file_handler)
            return params
        except IOError:
            logger.warn("'params.yaml' must be created before running a simulation.")

    params = get_config()
    name = getcwd().split('/')[-1]
    return name, params


def run(name, params):
    """ Spawn HDF file and create basic structure for the experiment.

    The `identifier` is composed from the name and a random hash value to
    differentiate different runs of the same simulation.
    """
    while True:             # ensure that the id is unique
        identifier = str(name) + '_' + hex(hash(random()))
        fname = identifier + '.hdf5'
        if not path.isfile(fname):
            break

    logger.info('\nSIMULATION STARTED >>> {}'.format(identifier))
    for k, v in params.items():
        logger.info('{}: {}'.format(k, v))

    with h5py.File(fname, 'x') as file_handler:
        mapping = eval('lambda x:' + params['function'])        # Care, eval is evil
        for value in params['values']:
            twinplane_set = mapping(value)
            twinplane_grp = file_handler.create_group('TWINS ' + str(twinplane_set))
            for sample in range(params['sample_size']):
                time = arrow.now().format('YYYY\'MM\'DD HH:mm:ss')
                logger.info('@' + time + ' Sampling... {}/{}'.format(
                    sample + 1, params['sample_size']))
                header, geometry = builder(twinplane_set, **params)
                sample_grp = twinplane_grp.create_group('sample{:04}'.format(sample))
                sample_grp.create_dataset('flake', data=geometry)
                if value == 0 and sample == 0:
                    file_handler.attrs['header'] = str(header)
    logger.info('SIMULATION ENDED >>> {}'.format(identifier))


def builder(tp, total_size=10000, sample_size=10, snapshot_intervall=1000,
        **kwargs):
    """ Generate a `Flake` instance from the `growth` module and sample the
    given parameters to generate meaningful statistics. The `geometry()` method
    of the flake is called at each snapshot, stored and returned with the corresponding
    header.
    """
    div, remainder = total_size // snapshot_intervall, total_size % snapshot_intervall
    thisFlake = Flake(*tp, **kwargs['flake'])
    header = thisFlake.geometry().keys()
    data = []

    logger.info('Twinplanes: {} \tTotal Size:{} \tSeed: {} \tTemp:[{}]'.format(
      thisFlake.twins, total_size, thisFlake.seed_shape, thisFlake.temp))

    for it in range(div):           # for all whole steps
        thisFlake.grow(snapshot_intervall)
        data.append(thisFlake.geometry().values())
    if remainder:                   # for the remaining ones
        thisFlake.grow(remainder)
        data.append(thisFlake.geometry().values())

    return (header, data)


def create():
    """ Create a new YAML file for the project.
    """


def main():
    parser = argparse.ArgumentParser(prog='FLagS')
    # parser.add_argument('command', choices=['run', 'create'])
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-c', '--create', action='store_const', const=True)
    group.add_argument('-r', '--run', action='store_const', const=True)
    # parser.add_argument('-r', '--run')
    args = parser.parse_args()

    if args.create is True:
        print('hooray')
    if args.run :
        print('trurun')
        # try:
            # os.environ['SUDO_USER']
            # install()
        # except KeyError:
            # sys.exit("Error: Must be run with sudo.")
    # elif args.command == 'aur':
        # loggy.warning("BEWARE still very experimental !")
        # try:
            # assert(os.environ['USER'] != 'root')
            # aur()
        # except AssertionError:
            # sys.exit("Error: Can not be run as superuser.\nExiting.")


#----------------
#  Boilerplate  -
#----------------
if __name__ == "__main__":
  main()
