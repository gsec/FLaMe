#!/usr/bin/env python2
## encoding: utf-8
# Statistics Collector
from __future__ import division
from os import path
from growth import Flake
from random import randint


def stats(iters=10**3, samples=10):
    TP_TRIALS = 5
    MODE = 'det'
    ID = randint(0, 9999)
    FNAME = 'stats_{}-[{}]_It:[{}K].txt'.format(ID, MODE, 0.001*iters)
    tp_func = lambda x: (-x, x)

    for r in range(TP_TRIALS):
        def collector():
            myFlake = Flake(*tp, seed='point', temp=0)
            myFlake.grow(iters, mode=MODE)
            myFlake.geometry()
            # content = '\n' + str(myFlake.geometry()) + '\n'
            # acc['string'].append(content)
            acc['aspect_ratio'] += myFlake.aspect_ratio
            acc['height'] += myFlake.height

        tp = tp_func(r)
        myF = Flake(*tp, seed='point', temp=0)
        acc = {'string': ['\n' + str(myF) + '\n' + 75*'='], 'aspect_ratio': 0,
                     'height': 0}
        folder = myF.daily_output()
        print("\nEvaluation with [{}] TPS!\t pUnique Ident -- {}\n".format(tp, ID))

        for s in range(samples):
            collector()
        mean_ar = acc['aspect_ratio'] / samples
        mean_h = acc['height'] / samples
        mean_msg = ('\n' + 50*'-' + "\nMEAN of [{}] -- \tAspect Ratio:"
                                "[{}]\tHeight:[{}]".format(samples, mean_ar, mean_h) +
                                "\n" + 50*'-' + "\n")
        with open(path.join(folder, FNAME), 'ab') as fil:
            fil.writelines(acc['string'])
            fil.write(mean_msg)


if __name__ == '__main__':
    stats(iters=10**5)
