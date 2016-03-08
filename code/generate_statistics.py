#!/usr/bin/env python2
## encoding: utf-8
# Statistics Collector
from __future__ import division
# from os import path
from growth import Flake
from random import randint
import arrow
import logging

## LOG CONFIG
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.FileHandler('statistics.log')
handler.setLevel(logging.INFO)
logger.addHandler(handler)


def stats(iters=10**3, samples=20):
    TP_TRIALS = 10
    MODE = 'det'
    ID = randint(0, 0xFFFF)
    tp_func = lambda x: (0, x)

    for r in range(TP_TRIALS):
        ar_acc, height_acc = 0.0, 0.0
        tp = tp_func(r)
        myF = Flake(*tp, seed='point', temp=0)
        time = '@' + arrow.now().format('YYYY\'MM\'DD HH:mm:ss')
        logger.info(str(myF) + '\n' + 80*'=')
        logger.info(time + " Sampling... Twinplanes: {} \tpUIdent:{}".format(tp, ID))

        for s in range(samples):
            myFlake = Flake(*tp, seed='point', temp=0)
            myFlake.grow(iters, mode=MODE)
            logger.info(myFlake.geometry())
            ar_acc += myFlake.aspect_ratio
            height_acc += myFlake.height

        mean_ar = ar_acc / samples
        mean_h = height_acc / samples
        mean_msg = ':: MEAN of [{}] :: \tAspect Ratio: {} \tHeight: {}\n'.format(
            samples, mean_ar, mean_h)
        logger.info(mean_msg)


if __name__ == '__main__':
    stats(iters=10**5)
