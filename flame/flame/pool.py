from flame.growth import Flake
import arrow
from multiprocessing import Pool
from random import random
import pylab

ROUNDS = 1000
FILENAME = 'tests/time_evals_{}.txt'.format(hex(hash(random())))

def flake_growth(flake, rounds=ROUNDS):
    flake.grow(rounds)

def control_group(runs):
    storm = [Flake(-x, x+1) for x in range(runs)]
    global start
    start = arrow.now()
    for each in storm:
        each.grow(ROUNDS)
    stop = arrow.now()
    return (stop-start).total_seconds()


def pool_tester(runs):
    storm = [Flake(-x, x+1) for x in range(runs)]
    start = arrow.now()
    with Pool() as p:
        p.map(flake_growth, storm)
    stop = arrow.now()
    return (stop-start).total_seconds()

if __name__ == '__main__':
    x = list(range(20))

    # control group
    cg_times = [control_group(arg) for arg in x]
    pylab.plot(x, cg_times)

    # in pool
    pg_times = [pool_tester(arg) for arg in x]
    pylab.plot(x, pg_times)

    mfg = """ Pool test log @ {time}

    Number of tasks:
{tasks}

    Control group:
{cg}

    Pool group:
{pg}
    """

    with open(FILENAME, 'w') as f:
        f.write(mfg.format(tasks=x, cg=cg_times, pg=pg_times, time=start))

    pylab.show()
