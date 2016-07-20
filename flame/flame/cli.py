"""
    FLaMe - a FLakeLAtticeMOdelEr
    **********************************

    This package is written in Python 2, nevertheless we seek maximum compatibility for
    Python 3. The main issue is the `mayavi` package not yet running properly in Python 3.
    Mayavi is required for plotting, although the simulations and all the other Flake()
    operations will also work under Python 3.
"""
from flame import simulation as sim
import os, argparse


def main():
    parser = argparse.ArgumentParser(prog='FLaMe')
    parser.add_argument('command', choices=['run', 'create'])
    parser.add_argument("name", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    if args.command == 'run':
        sim.run(*sim.read_params())
    if args.command == 'create':
        if not args.name:
            args.name = [os.getcwd().split('/')[-1]]
        sim.create('_'.join(args.name))
