"""
    FLaMe - a FLakeLAtticeMOdelEr
    **********************************
    - command line entry point -

    This package is written in Python 2, nevertheless we seek maximum compatibility for
    Python 3. The main issue is the `mayavi` package not yet running properly in Python
    3. Mayavi is required for 3D single flake plotting, although the simulations and all
    the other Flake() operations will also work under Python 3.
"""
from flame import simulation, paint
import os, argparse


def main():
    parser = argparse.ArgumentParser(prog='FLaMe')
    parser.add_argument('command', choices=['run', 'create', 'paint'])
    parser.add_argument("name", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    if args.command == 'run':
        simulation.run(simulation.get_params())
    if args.command == 'create':
        if not args.name:
            args.name = [os.getcwd().split('/')[-1]]
        simulation.create('_'.join(args.name))
    if args.command == 'paint':
        paint.mean_plot(*args.name)
