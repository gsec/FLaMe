#!/usr/bin/env python
# encoding: utf-8

from os import path
from growth import Flake
from random import random


def sample(tp=tuple(), iter=2*10**5, noise=1e-5, runs=10, export=False):
  ident = str(hash(random()))[-4:]
  folder = Flake().daily_output()
  fname = 'tRun_It:[{}K]_Noise:[{}]_Rounds:[{}]_{}.txt'.format(
    iter/1000, noise, runs, ident)
  print("twiNs:{}_".format(tp) + fname[:-4] + 'id')

  for r in range(runs):
    print("\nRUN:\t{}\n{}".format(r, 70*'.'))
    obj = Flake()
    obj.grow(iter, noise)
    content = '\n' + str(obj.geometry()) + '\n'

    with open(path.join(folder, fname), 'ab') as fil:
      if r == 0:
        fil.write("{} -- NOIZ:\t[{}]{}{}{}".format(
          obj, noise, '\n', 75*'=', '\n'))
      fil.write(content)
  if export:
    obj.export('COORDS_' + fname[:-4])

if __name__ == '__main__':
  sample(tp=(-5, 3), iter=10**8, runs=100)
