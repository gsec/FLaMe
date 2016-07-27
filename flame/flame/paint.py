# coding: utf-8

from os import path, mkdir
import logging
import pandas as pd
from bokeh.plotting import figure, output_file, show
from flame.settings import GRAPH_OUTPUT

logger = logging.getLogger(__name__)


def mscatter(p, x, y, marker, color):
    p.scatter(x, y, marker=marker, size=5,
              line_color="navy", fill_color=color, alpha=0.5)


def mean_plot(fname):

    print(fname)
    name = fname.split('.')[0]
    fullpath = path.join(GRAPH_OUTPUT, name)

    try:
        mkdir(fullpath)
    except OSError as e:
        if e.errno == 17:
            logger.debug('Did not create {}.\n{}'.format(fullpath, e))

    p = figure(title=name)
    p.grid.grid_line_color = "white"
    p.background_fill_color = "#eeeeee"

    chosen_col = 'height'

    with pd.HDFStore(fname, 'r') as h5:
        twins = [group for group in h5.root if 'twin' in str(group)]
        for idxt, twin in enumerate(twins):
            tmp = []
            for sample in twin:
                loc = str(sample).split()[0]
                flake = h5.select(loc)
                tmp.append(flake)
            flake_sum = pd.concat(tmp)
            group_by_index = flake_sum.groupby(flake_sum.index)
            mean_flake = group_by_index.mean()

            attrib = mean_flake[chosen_col]
            iters = mean_flake['iter']


            norm = len(twins)-1
            if norm == 0:
                logger.warn("Single Flakes, no averaging.")
                norm = 1
            red = 255*idxt/norm
            green = 0
            blue = 255*(1-idxt/norm)
            logger.info("{} RGB: {} {} {}".format(str(twin), red, green, blue))

            mscatter(p, iters, attrib, "circle", (red, green, blue))
    output_file(path.join(fullpath, chosen_col + '.html'))
    show(p)
