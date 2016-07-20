# from bokeh.plotting import figure, output_file, show
import tables as tb
import logging
from flame import simulation as sim

logger = logging.getLogger(__name__)


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
    try:
        import arrow
    except ImportError:
        return "n0timet0day"
    return arrow.now().format('YYYY\'MM\'DD HH:mm:ss')


params = {
    'flake': {'temp': 40, 'trail_length': 20},
    'function': 'set(2*y for y in range(-x+1, x))',
    'name': 'clori',
    'sample_size': 5,
    'snapshot_intervall': 1000,
    'total_size': 200000,
    'values': [1, 2, 0, 3]}

def write_data(fname):

    mapping = eval('lambda x:' + params['function'])        # Care, eval is evil

    with tb.open_file(fname, mode='w', title='Test file') as h5file:
        for each_value in params['values']:
            twinplane_set = mapping(each_value)
            for sample in range(params['sample_size']):
                logger.info(' @{time} Twinplanes {twins} Total Size: {size}'.format(
                    time=get_time(), twins=twinplane_set, size=params['total_size']))
                flake_gen = sim.builder(twinplane_set, **params)
                table = h5file.create_table(
                    '/twins_{}'.format('_'.join([str(tp) for tp in twinplane_set])),
                    'flake{:03}'.format(sample), Stage, 'FlakeSample',
                    createparents=True)
                snapshot = table.row
                for stage in flake_gen:
                    for k, v in stage.items():
                        snapshot[k] = v
                    snapshot.append()
                table.flush()


logger.setLevel(logging.DEBUG)
write_data('testname')
