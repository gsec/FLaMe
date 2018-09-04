""" Time benchmarking script.

Compares runtime between single-threaded and multi-threaded Flake growth simulation.
"""
import arrow
import pylab
import logging
from os import path
from numpy import mean
from random import random
from multiprocessing import Pool

from flame.growth import Flake
from flame.settings import TIMING_OUTPUT

logger = logging.getLogger()


def flake_growth(flake, itersize=1000):
    """ This is a wrapper method which we can hand in Flake instances. Each of those
    instances will be processed parallely by `multithread.Pool`.
    """
    flake.grow(itersize)


def snow_storm(flake_count, parallel=False, **kwargs):
    """ First create a whole bunch of Flakes (one might say a snow storm seeded in a
    cloud), which are then passed to the `flake_growth()` wrapper. If we are in parallel
    mode this gets executed on multiple threads by starmap(). The twinplanes are a simple
    function, slowly varying, to reduce effects of a special case.
    """
    cloud = [(Flake(-n, n+1), kwargs['itersize']) for n in range(flake_count)]
    logger.info("Generate cloud of size [{}]: parallel={}".format(flake_count, parallel))
    start = arrow.now()

    if parallel:
        with Pool() as p:
            p.starmap(flake_growth, cloud)
    else:
        for flake_seed in cloud:
            flake_growth(*flake_seed)

    eval_time = (arrow.now() - start).total_seconds()
    logger.info("Took total of {} seconds.\n".format(eval_time))
    return eval_time


def main(**kwargs):
    """ Here we catch missing arguments with very conservative numbers. This is only
    intended for testing etc. the correct parameters should be passed in the call.
    """
    try:
        itersize = kwargs['itersize']
    except KeyError:
        kwargs['itersize'] = itersize = 500
    try:
        benchmark_range = kwargs['benchmark_range']
    except KeyError:
        kwargs['benchmark_range'] = benchmark_range = 5

    benchmarks = list(range(1, benchmark_range + 1))
    logging.info("Benchmark Range: {}\n".format(benchmark_range))

    FILENAME = 'time_evals_#{}_@{}it_{}.txt'.format(benchmark_range, itersize,
                                                    hex(hash(random())))
    outpath = path.join(TIMING_OUTPUT, FILENAME)

    # control group
    cg_times = [snow_storm(those, parallel=False, **kwargs) for those in benchmarks]
    pylab.plot(benchmarks, cg_times)

    # in pool
    pg_times = [snow_storm(those, parallel=True, **kwargs) for those in benchmarks]
    pylab.plot(benchmarks, pg_times)

    message_skel = """Pool test log @ {time}
We build an incrementing batch of flakes and let them grow.  The *control group* is
called normally, i.e. with a single thread, the *pool group* is called by the starmap()
function of multithreading.Pool().

    Control group timings:
{cg}

    Pool group timings:
{pg}

    Timing ratios (control/pool):
{cpg}

Total number of flakes in last cloud: {tasks}
Iterations per flake growth: {it}

Total speedup (averaged time ratios):
          -------------
        ==  {avg:.3f}  ==
          -------------
"""
    t_ratios = [(cg/pg) for (cg, pg) in zip(cg_times, pg_times)]
    t_avg = mean(t_ratios)

    with open(outpath, '+x') as f:
        f.write(message_skel.format(tasks=benchmark_range, cg=cg_times, pg=pg_times,
                                    time=arrow.now(), cpg=t_ratios, avg=t_avg,
                                    it=itersize))
    pylab.show()


if __name__ == '__main__':
    SIZE = 100000
    BENCHMARKS = 20

    main(itersize=SIZE, benchmark_range=BENCHMARKS)
