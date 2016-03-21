#!/usr/bin/env python2
## encoding: utf-8
# Statistics Collector
from __future__ import division
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
    tp_func = lambda x: (-x, x + 1)
    logger.info(39*'.-' + '.')
    logger.info('Flake Statistics -- Mode: [{}] => {} samples á {}K atoms\n'.format(
       MODE.upper(), samples, iters/1000.))

    for r in range(TP_TRIALS):
        ar_acc, height_acc = 0.0, 0.0
        tp = tp_func(r)
        myF = Flake(*tp, seed='point', temp=0)
        time = '@' + arrow.now().format('YYYY\'MM\'DD HH:mm:ss')
        logger.info(time + ' Sampling... ')
        logger.info('Twinplanes: {} \tpUIdent:{} \tSeed: {} \tTemp:[{}]'.format(
          tp, ID, myF.seed_shape, myF.temp))

        for s in range(samples):
            myFlake = Flake(*tp, seed='point', temp=0)
            myFlake.grow(iters, mode=MODE)
            logger.info('Aspect Ratio: {aspect_ratio:>5.2f}, Area:{area:>9.2f}, '
                        'Radius: {radius:>5.2f}, Height: {height:>5.2f}'.format(
                          **myFlake.geometry()))
            ar_acc += myFlake.aspect_ratio
            height_acc += myFlake.height

        mean_ar = ar_acc / samples
        mean_h = height_acc / samples
        mean_msg = (':: MEAN of [{}] :: \tASPECT RATIO: {:>5.2f} \tHEIGHT: '
                    '{:>5.2f}\n').format(samples, mean_ar, mean_h)
        logger.info(mean_msg)


if __name__ == '__main__':
    stats(iters=10**6, samples=100)
