"""
    FLaMe - a FLakeLAtticeMOdelEr
    **********************************
    - command line entry point -

    A hcc/fcc crystal growth simulation, supporting different growth modes, arbitrary
    twinplanes and an interface for bulk simulation and statistics.

    Support of Python2 is no longer a goal, we transfer it to Py3 while dropping mayavi
    and instead exporting the atom coordinates in a proper blender format.
"""
from flame import simulation, paint
import os, argparse


def main():
    parser = argparse.ArgumentParser(prog='FLaMe')
    parser.add_argument('command', choices=['run', 'create', 'paint', 'extract'])
    parser.add_argument("name", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    if args.command == 'run':
        simulation.run(simulation.get_params())
    elif args.command == 'create':
        if not args.name:
            args.name = [os.getcwd().split('/')[-1]]
        simulation.create('_'.join(args.name))
    elif args.command == 'paint':
        paint.mean_plot(*args.name)
    elif args.command == 'extract':
        for fname in args.name:
            paint.extractor(fname)
