#!/usr/bin/env python2
## encoding: utf-8
# Statistics Collector
from __future__ import division
from growth import Flake
# from random import randint
# from mayavi import mlab
import arrow
import logging

## LOG CONFIG
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.FileHandler('movies.log')
handler.setLevel(logging.INFO)
logger.addHandler(handler)


# @mlab.show
# @mlab.animate(ui=False)
# def animate(sampling=20):
#     f = Flake(-2, 2, trail=sampling)
#     coords  = f.plot(pipeline=True)
#
#     src = mlab.points3d(*coords, colormap="gist_ncar", scale_factor=0.1, vmin=0, vmax=15)
#     ms = src.mlab_source
#     fig = mlab.gcf()
#
#     while True:
#         f.grow(sampling, mode='det')
#         f.carve()
#         x, y, z, c = f.plot(pipeline=True)
#         ms.reset(x=x, y=y, z=z, scalars=c)
#         mlab.yaw(10)
#         fig.scene.reset_zoom()
#         yield
#
# a = animate()


# def anim():
    # """Animate."""
    # while 1:
        # f.grow(TL)
        # f.carve()
        # x, y, z, c = f.plot(pipeline=True)
        # ms.reset(x=x, y=y, z=z, scalars=c)
        # mlab.yaw(10)
        # fig.scene.reset_zoom()
        # yield


# ----------
# -  main  -
# ----------
def animate(sampling=20):
    f = Flake(-2, 2, trail=sampling)
  # f = Flake(size=71)
  # atomic_num = (2*f.seed_size + 1)**3
  # f.tag = tag

    all_timings = '\n'
  # f.plot(save=True, tag='seed')
    for round in range(10):
        # atomic_num += binning
        g_start = arrow.now()
        f.grow(sampling)
        g_end = arrow.now()
        p_start = arrow.now()
        f.plot(save=True)
        p_end  = arrow.now()
        p_delta = (p_end - p_start).total_seconds()
        g_delta = (g_end - g_start).total_seconds()
        timing_string = ("Count: {} atoms. Added {} atoms in {} "
                        "sec and rendered in {} sec\n").format(
                            f.iter, sampling, g_delta, p_delta)
        # qprint(timing_string, quiet=Q)
        all_timings += timing_string
    # fname = path.join(f.daily_output(), 'timings_' + '.txt')
    # with open(fname, 'a+') as tfile:
        # tfile.write(f.__repr__() + all_timings + '\n')

if __name__ == '__main__':
    animate()

 # def main():
   # for x in range(10):
     # tag = 'pro001_' + str(x)
     # binn = (1 + x) * 50
     # animate(tag, binning=binn)

#  if __name__ == '__main__':
#    main()

# def stats(iters=10**3, samples=20):
    # TP_TRIALS = 10
    # MODE = 'det'
    # ID = randint(0, 0xFFFF)
    # tp_func = lambda x: (-x, x + 1)
    # logger.info(39*'.-' + '.')
    # logger.info('Flake Statistics -- Mode: [{}] => {} samples á {}K atoms\n'.format(
       # MODE.upper(), samples, iters/1000.))

    # for r in range(TP_TRIALS):
        # ar_acc, height_acc = 0.0, 0.0
        # tp = tp_func(r)
        # myF = Flake(*tp, seed='point', temp=0)
        # time = '@' + arrow.now().format('YYYY\'MM\'DD HH:mm:ss')
        # logger.info(time + ' Sampling... ')
        # logger.info('Twinplanes: {} \tpUIdent:{} \tSeed: {} \tTemp:[{}]'.format(
          # tp, ID, myF.seed_shape, myF.temp))

        # for s in range(samples):
            # myFlake = Flake(*tp, seed='point', temp=0)
            # myFlake.grow(iters, mode=MODE)
            # logger.info('Aspect Ratio: {aspect_ratio:>5.2f}, Area:{area:>9.2f}, '
                        # 'Radius: {radius:>5.2f}, Height: {height:>5.2f}'.format(
                          # **myFlake.geometry()))
            # ar_acc += myFlake.aspect_ratio
            # height_acc += myFlake.height

        # mean_ar = ar_acc / samples
        # mean_h = height_acc / samples
        # mean_msg = (':: MEAN of [{}] :: \tASPECT RATIO: {:>5.2f} \tHEIGHT: '
                    # '{:>5.2f}\n').format(samples, mean_ar, mean_h)
        # logger.info(mean_msg)


# if __name__ == '__main__':
    # stats(iters=10**6, samples=100)
