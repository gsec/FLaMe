# coding: utf-8

from os import path, mkdir
import logging
import pandas as pd
from bokeh.plotting import figure, output_file, show
from flame.settings import GRAPH_OUTPUT, HDF_METADATA, get_colors

logger = logging.getLogger(__name__)


def mscatter(p, x, y, marker, color, legend=None):
    """ Shamelessly stolen from bokeh documentation.
    """
    p.scatter(x, y, marker=marker, size=5,
              line_color=color, fill_color=color, alpha=0.6, legend=legend)


def mean_plot(file_path, *columns, **kwargs):
    """ Paint bokeh over the mean of samples for each twinplane configuration.

    Outputs to html file in `output/graph/<filebasename>/<column>.html`.
    """
    fname = path.basename(file_path)
    name = fname.split('.')[0]
    output_dir = path.join(GRAPH_OUTPUT, name)
    try:
        mkdir(output_dir)
    except OSError as e:
        if e.errno == 17:
            logger.debug('Did not create {}.\n{}'.format(output_dir, e))

    if not columns:
        columns = ('height',)

    for col in columns:
        figtitle = "{} ({})".format(name, col)
        p = figure(title=figtitle, width=1024, height=768)
        p.grid.grid_line_color = "white"
        p.background_fill_color = "#eeeeee"

        for plane, twin, color, mean_flake in averaged_planes(fname):
            Y = mean_flake[col]
            X = mean_flake['iter']
            mscatter(p, X, Y, "circle", color, legend=plane)

        p.legend.location = 'bottom_right'
        output_file(path.join(output_dir, col + '.html'))
        show(p)

"""
TODO:
    separate averaging from plotting

    create buttons for each twin plane as selector to switch their display

    each tp selector should also have a switch to display all or only averaged results
"""


def averaged_planes(fname, choices=None):
    with pd.HDFStore(fname, 'r') as h5:
        twins  = [group for group in h5.root if HDF_METADATA not in str(group)]
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
            try:
                plane = str(h5.select('/parameters')['planes'][index]).replace(
                    'set', 'Twinplanes: ')
            except KeyError:
                plane = str(twin).split()[0]
            yield (plane, twin, color, mean_flake)
