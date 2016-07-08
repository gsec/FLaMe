#!/usr/bin/env python2
# encoding: utf-8
"""                           FLagS - a FLakeLatticeGrowthSimulator

                              <guilherme.stein@physik.uni-wuerzburg.de>

Define the general Structure of the Simulations (here we use this as a synonym to
experiments). We will create a hierarchy in the HDF data format and group those elements
structure for the experiments we want to run with our basic flake growth module
"""

from __future__ import print_function, division, generators
from flame.growth import Flake
from random import random
from os import getcwd, path
import tables as tb
import logging
import yaml
import sys

# LOG CONFIG
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# canonical file name for configs
PARAMS_YAML = 'sim_params.yaml'


class Stage(tb.IsDescription):
    """ Definition class for output of `geometry()` at a snapshot.
    """
    radius          = tb.Float32Col()
    height          = tb.Float32Col()
    aspect_ratio    = tb.Float32Col()
    area            = tb.Float32Col()
    layers          = tb.Int32Col()
    iter            = tb.Int32Col()


def get_time():
    """ Returns formatted time string.

    If arrow import fails, warning is issued, and dummy string is returned.
    """
    try:
        import arrow
    except ImportError:
        logger.warn("Import of arrow failed. No time is recorded. Please check "
                    "dependencies.")
        return "n0timet0day"
    return arrow.now().format('YYYY\'MM\'DD HH:mm:ss')


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

    logger.info('\t\t STARTED >>> {}'.format(identifier))
    for k, v in params.items():
        logger.info('{}: {}'.format(k, v))

    mapping = eval('lambda x:' + params['function'])        # Care, eval is evil

    with tb.open_file(fname, mode='w', title='Test file') as h5file:
        for each_value in params['values']:
            twinplane_set = tuple(set(mapping(each_value)))
            twinplane_string = '/twins_{}'.format(
                '_'.join('0' + str(abs(tp)) if tp < 0 else str(tp) for tp in
                         sorted(list(twinplane_set))))
            logger.info(' @{time} Twinplanes {twins} Total Size: {size}'.format(
                time=get_time(), twins=twinplane_set, size=params['total_size']))

            for sample in range(params['sample_size']):
                table = h5file.create_table(
                    twinplane_string, 'flake{:03}'.format(sample), Stage, 'FlakeSample',
                    createparents=True)
                logger.info('Sampling... {}/{}'.format(
                    sample + 1, params['sample_size']))

                snapshot = table.row
                flake_gen = builder(twinplane_set, **params)

                for stage in flake_gen:
                    for k, v in stage.items():
                        snapshot[k] = v
                    snapshot.append()
                table.flush()
    logger.info('SIMULATION ENDED >>> {}'.format(identifier))


def store_hashed_params():
    #TODO: put attributes into hdf, ensure sim_params.yaml is not altered afterwards.
    pass

def read_params():
    """ Ensures the parameters are read in from a YAML file.

    A simulation is defined through its configuration file. Here we use YAML
    since it is human- and machine-readable.
    """
    def get_config():
        try:
            with open(PARAMS_YAML, 'r') as file_handler:
                params = yaml.load(file_handler)
            return params
        except IOError:
            sys.exit("'{}' must be created before running a simulation. "
                     "Try the `create` command.".format(PARAMS_YAML))

    params = get_config()
    name = getcwd().split('/')[-1]
    return name, params


def builder(tp, total_size=10000, snapshot_intervall=1000, **kwargs):
    """ Generate a `Flake` instance from the `growth` module and sample the
    given parameters to generate meaningful statistics. The `geometry()` method
    of the flake is called at each snapshot, stored and returned with the corresponding
    header.
    """
    try:
        xargs = kwargs['flake']
    except KeyError:
        xargs = {'': None}
    thisFlake = Flake(*tp, **xargs)

    while thisFlake.iter < total_size:
        thisFlake.grow(snapshot_intervall)
        yield thisFlake.geometry()


def create(project_name):
    """ Create a default YAML file for the project.

    Optional argument: @name   :: If not provided, naming defaults to directory name.
    """
    if path.isfile(PARAMS_YAML):
        sys.exit("Parameters file already exists in {}\n".format(getcwd()) +
                 "Please create each simulation in a separate folder.")

    skeleton = """\
# Configuration file generated by FLaMe @{time}.
# This file should be in its own directory where the simulation is started.
# Parameters specified here will be used for the growth simulation. Feel free
# to edit responsibly.

    # Name of the Project
name: {name}
    # Function string will be parsed as 'lambda x: ' + function
function: "(0, x)"
    # Values that will be fed as input to f(x)
values: [0, 1, 2, 5]
    # How many averaging runs are done with exact same parameters
sample_size: 5
    # Interval between Flake.geometry() calls
snapshot_intervall: 100
    # End size of the flake in atoms
total_size: 100000

    # These options are passed to growth.Flake() instances.
flake:
    temp: 50
    seed: 'point'
    trail: 20\
"""
    output = skeleton.format(name=project_name, time=get_time())
    with open(PARAMS_YAML, 'w') as handler:
        handler.write(output)
    logger.info('Parameter file written to {}'.format(PARAMS_YAML))
