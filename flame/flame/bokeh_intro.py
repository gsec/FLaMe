# coding: utf-8

from bokeh.charts import Line, output_file, show
from flame.settings import GRAPH_OUTPUT
from os import path, mkdir
import tables as tb


def simple_plot(fname):

    print(fname)
    name = fname.split('.')[0]
    fullpath = path.join(GRAPH_OUTPUT, name)

    try:
        mkdir(fullpath)
    except OSError:
        print('couldnt make IITTT')

    with tb.open_file(fname, 'r') as f:
        single = f.root.twins_04_02_0_2.flake002

        output_file(path.join(fullpath, 'aspects' + '.html'))

        aspect = single.col("aspect_ratio")

        show(Line(aspect))
