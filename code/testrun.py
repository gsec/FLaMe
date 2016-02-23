from sys import exit
from os import path
from growth import Flake
from random import random


def testflake(tp=tuple(), iterate=10**5, noise=1e-4, runs=10, export=False):
  ident = str(hash(random()))[-4:]

  for r in range(runs):
    print("RUN:\t{}".format(r))
    obj = Flake()
    obj.grow()
    obj.grow(iterate, noise)

    content = '\n' + str(obj.geometry()) + '\n'
    folder = obj.daily_output()
    fname = 'tRun_It:[{}K]_Noise:[{}]_Rounds:[{}]-{}.txt'.format(
      iterate/1000, noise, runs, ident)

    with open(path.join(folder, fname), 'ab') as fil:
      if r == 0:
        # fil.write(str(obj) + '\n' + 75*'=' + '\n')
        fil.write("{} -- NOIZ:\t[{}]{}{}{}".format(obj, noise, '\n', 75*'=', '\n'))  # + + + '\n')
      fil.write(content)
  if export:
    obj.export('COORDS_' + fname[:-4])


def main(*args, **kwargs):
  try:
    return testflake(*args, **kwargs)
  except KeyboardInterrupt:
    print("TestRun has been cancelled by user!")
    exit(0)
