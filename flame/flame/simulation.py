import sys
import logging
import pandas as pd
from random import random
from os import getcwd, path
from multiprocessing import Pool

from flame.growth import Flake
from flame.settings import HDF_EXT, PARAMS_YAML, get_time, get_skel, get_params

logger = logging.getLogger(__name__)


def create(project_name):
    """ Create a default YAML file for the project.

    Optional argument: @name   :: If not provided, naming defaults to directory name.
    """
    if path.isfile(PARAMS_YAML):
        sys.exit("Parameters file already exists in {}\n".format(getcwd()) +
                 "Please create each simulation in a separate folder.")

    output = get_skel(project_name)

    with open(PARAMS_YAML, 'w') as handler:
        handler.write(output)
    logger.info('Parameter file written to {}'.format(PARAMS_YAML))


def tp_gen(params):
    tp_func = eval('lambda x:' + params['function'])        # Care, eval is evil

    if type(params['values']) is str:
        objct = eval(params['values'])
    else:
        objct = params['values']

    mapping = map(tp_func, objct)
    return [set(x) for x in mapping]


def growth_sample(twin, params):
    return pd.DataFrame(builder(twin, **params))


def run(params=None):
    """ Spawn HDF file and create basic structure for the experiment.

    The `identifier` is composed from the name and a random hash value to
    differentiate different runs of the same simulation.

    * create an unique identifier value
    * create a mapping through lambda evaluation
    * generate twin plane configuration for each value

    New runners function:

    * with pandas
    * meta attributes
    * proper loops
    """
    if not params:
        params = get_params()
    identifier = params['name'] + '_' + hex(hash(random()))
    fname = identifier + HDF_EXT

    logger.info('STARTED >>> {} @ {}'.format(identifier, ' :: '.join(get_time())))
    for k, v in params.items():
        logger.info('\t\t{}: {}'.format(k, v))

    params['twins'] = twins = tp_gen(params)

    with pd.HDFStore(fname, title=identifier) as h5:
        h5.put('/parameters', pd.Series(str(params)))

        for tp_idx, twin in enumerate(twins):
            twin_loc = 'twinplane{:02}'.format(tp_idx)
            logger.info('F>> TP: {tp}  Samples:{sm}  @{tm}  Total Size: {sz}'.format(
                tp=twin,
                sm=params['sample_size'],
                tm=get_time()[1],
                sz=params['total_size']))

            # here comes the data crunching
            args = ((twin, params) for _ in range(params['sample_size']))
            with Pool() as p:
                samples = p.starmap(growth_sample, args)

            # and the storage on to disk
            for idx, sample in enumerate(samples):
                flake_loc = 'flake{:03}'.format(idx)
                location = "/".join((twin_loc, flake_loc))
                h5.put(location, sample)

    logger.info('ENDED >>> {} @ {}'.format(identifier, ' :: '.join(get_time())))


def builder(tp, total_size=10000, snapshot_interval=1000, **kwargs):
    """ Create generator that yields the geometry of growing Flake.

    Generate a `Flake` instance from the `growth` module. By consuming an item we grow
    the flake by the amount in `snapshot_interval` and yield the output of our
    `geometry` method.
    """
    try:
        xargs = kwargs['flake']
    except KeyError:
        xargs = {'': None}

    thisFlake = Flake(*tp, **xargs)

    while thisFlake.iter < total_size:
        thisFlake.grow(snapshot_interval)
        yield thisFlake.geometry()
    if 'name' in kwargs:
        thisFlake.carve()
        thisFlake.export_coordinates(kwargs['name'])
