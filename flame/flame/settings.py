""" FLaMe settings
"""
import sys
import yaml
from os import environ, path
from matplotlib.cm import viridis
from logging import getLogger


logger = getLogger(__name__)

""" This file stores global settings used across all simulations.

The options here are set as how the program behaves, while in `sim_params.yaml` we define
simulations specific parameters.

Output: It is easiest to set an environment variable and fLame outputs data relative to
that path. `$ export FLAME_OUTPUT=/your/path/here`

Distance => `DIFF_CAP`: This variable denotes the maximum distance (in atomic radii)
until which atoms are considered nearest neighbors when checked. For touching atoms this
is two atomic radii, while we give it here 2.1 to include small deviations due to
numerical rounding errors when calculating the distance.
"""

# These are the variables we are setting in this file and are exporting to the simulation.
global FLAME_OUTPUT, GROW_OUTPUT, GRAPH_OUTPUT
global DIFF_CAP, PICKLE_EXT, HDF_EXT

try:
    FLAME_OUTPUT = environ['FLAME_OUTPUT']
except KeyError as e:
    logger.warn("No {} specified. Path is chosen relative to `settings.py`.".format(e))
    fallback = path.dirname(path.dirname(path.abspath(__file__)))
    FLAME_OUTPUT = path.join(fallback, 'output')
logger.info("\nFLaMe output path set to: {}".format(FLAME_OUTPUT))

GROW_OUTPUT = path.join(FLAME_OUTPUT, 'grow')
SIM_OUTPUT = path.join(FLAME_OUTPUT, 'sim')
GRAPH_OUTPUT = path.join(FLAME_OUTPUT, 'graph')
DIFF_CAP = 2.1      # This is the maximum distance between two atoms to be considered
                    # nearest neighbors.
PICKLE_EXT = '.flm'
HDF_EXT = '.h5'
HDF_METADATA = 'parameters'
PARAMS_YAML = 'sim_params.yaml'


def get_time():
    """ Returns a list with two items: date and time.
    Both are formatted as unicode strings.

    If arrow import fails, warning is issued, and dummy string is returned.
    """
    try:
        import arrow
    except ImportError:
        logger.warn("Import of arrow failed. No time is recorded. Please check "
                    "dependencies.")
        return "n0timet0day"
    return arrow.now().format('YYYY-MM-DD HH-mm-ss').split(' ')


def get_skel(project_name):
    """ Returns formatted default yaml file.
    """

    skeleton = """\
# Configuration file generated by FLaMe @{time}.
# This file should be in its own directory where the simulation is started.
# Parameters specified here will be used for the growth simulation. Feel free
# to edit responsibly.

    # Name of the Project
name: '{name}'
    # Function string will be parsed as 'lambda x: ' + function
function: '{TP_FUNC}'
    # Values that will be fed as input to f(x)
values: {VALS}
    # How many averaging runs are done with exact same parameters
sample_size: {SMPSIZE}
    # Interval between Flake.geometry() calls
snapshot_interval: {SNAPSHOT}
    # End size of the flake in atoms
total_size: {TOTAL_SIZE}

    # These options are passed to growth.Flake() instances.
flake:
    temp: {FLAKE_TEMP}
    seed: '{FLAKE_SEED}'\
"""
    # defining defaults for sim_params file
    rendered = skeleton.format(
                    time=get_time(),
                    name=project_name,
                    TP_FUNC="(-(x%2), x, -x**2)",
                    VALS=[0, 1, 2, 3, 5],
                    SMPSIZE=10,
                    SNAPSHOT=100,
                    TOTAL_SIZE=200000,
                    FLAKE_TEMP=30,
                    FLAKE_SEED='point')
    return rendered


def get_params():
    """ Ensures the parameters are read in from a YAML file.

    A simulation is defined through its configuration file. Here we use YAML
    since it is human- and machine-readable.
    """
    try:
        with open(PARAMS_YAML, 'r') as file_handler:
            params = yaml.load(file_handler)
        return params
    except IOError:
        sys.exit("'{}' must be created before running a simulation. "
                    "Try the `create` command.".format(PARAMS_YAML))


def get_colors(twins):
    """ Generate colors for different twinplanes.

    Scales the index evenly between twinplanes and returns a RGBA float value from the
    colormap (here: viridis).
    """
    norm = float(len(twins) - 1)
    if norm == 0:
        logger.warn("Single Flakes, no averaging.")
        norm = 1

    color_filler = []
    for idx, tp in enumerate(twins):
        scala = idx/norm
        rescaled = tuple(x*255 for x in viridis(scala)[:-1])
        color_filler.append(rescaled)

    return zip(range(len(twins)), twins, color_filler)
