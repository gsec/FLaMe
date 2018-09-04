import logging
import pandas as pd
from os import path, mkdir
from bokeh.layouts import column as bokeh_column
from bokeh.plotting import figure, output_file, show

from flame.growth import Flake
from flame.settings import GRAPH_OUTPUT, HDF_METADATA, get_colors

logger = logging.getLogger(__name__)


def extractor(file_path):
    h5 = pd.HDFStore(file_path, 'r')

    logger.info("{}\t{}\t{}".format(20*'>', file_path, 20*'<'))
    for key, val in h5.parameters.items():
        logger.info("{}:{}{}".format(key, (20-len(key))*" ", val))
    h5.close()


def mscatter(p, x, y, marker, color, legend=None):
    """ Shamelessly stolen from bokeh documentation.

    This is basically an intermediate layer to simplify scatter creation.
    """
    p.scatter(x, y, marker=marker, size=5,
              line_color=color, fill_color=color, alpha=0.6, legend=legend)


def cols_output(file_path):
    """ Create path and return columns
    """
    fname = path.basename(file_path)
    name = fname.split('.')[0]
    output_dir = path.join(GRAPH_OUTPUT, name)
    try:
        mkdir(output_dir)
    except OSError as e:     # pragma: nocover
        if e.errno == 17:
            logger.info('Directory exists. Did not create {}.\n{}'.format(output_dir, e))
        elif e.errno == 2:
            logger.warning('Invalid path! Did not create {}.\n{}'.format(output_dir, e))
        else:
            raise e
    return output_dir, fname, name


def mean_plot(file_path, *columns):
    """ Paint bokeh over the mean of samples for each twinplane configuration.

    Outputs to html file in `output/graph/<filebasename>/geometry_##_cols.html`.
    """
    output_dir, fname, name = cols_output(file_path)

    if not columns:
        columns = ('aspect_ratio',)
    elif 'all' in columns:
        columns = list(Flake().geometry().keys())
        columns.pop(-2)         # remove 'iter' variable from 'all'

    logger.info("Generating columns:\t{}".format(columns))
    plot_cols = []

    for col in columns:
        figtitle = "{} ({})".format(name, col)
        p = figure(title=figtitle, width=1024, height=768)
        p.grid.grid_line_color = "white"
        p.background_fill_color = "#eeeeee"

        for plane, color, mean_flake in averaged_planes(fname):
            Y = mean_flake[col]
            X = mean_flake['iter']
            mscatter(p, X, Y, "circle", color, legend=plane)

        p.legend.location = 'bottom_right'
        plot_cols.append(p)

    out_fname = 'geometry_{}_cols.html'.format(len(columns))
    output_file(path.join(output_dir, out_fname))
    show(bokeh_column(plot_cols))


def averaged_planes(fname, choices=None):
    """
    Extracts data from HDF file and averages over the flakes with same configuration.

    TODO:
        * separate averaging from plotting
        * create buttons for each twin plane as selector to switch their display
        * each tp selector should also have a switch to display all or only averaged
          results
    """
    with pd.HDFStore(fname, 'r') as h5:
        twins = [group for group in h5.root if HDF_METADATA not in str(group)]
        if choices:
            twins = [x for i, x in enumerate(twins) if i in choices]
        for index, twin, color in get_colors(twins):
            tmp = []
            for sample in twin:
                loc = str(sample).split()[0]
                flake = h5.select(loc)
                tmp.append(flake)

            flake_sum = pd.concat(tmp)
            group_by_index = flake_sum.groupby(flake_sum.index)
            mean_flake = group_by_index.mean()
            plane = str(h5.select('/' + HDF_METADATA)['twins'][index]).replace(
                'set', 'Twinplanes: ')
            yield (plane, color, mean_flake)
