import os, argparse
from flame import simulation, paint


def main():
    """ Command line entry point.
    """
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
