import unittest
from flame.paint import get_planes
from os import path


# First start with HDF file access

# small simulation as test file
p = path.dirname(path.abspath(__file__))


class TestAccess(unittest.TestCase):
    # if we ran a simulation, show what twinplanes are in it
    def test_listing(self):
        try:
            h5 = path.join(p, 'sim_test.h5')
            listing = get_planes(h5)
            print(listing)
        except:
            print("No testfile found! No tests performed in {}!\n".format(__file__))


# We want to be able to select a twinplane by name (tuple, index?)
# tp = select_plane(h5, (-1, 1))

if __name__ == '__main__':
    unittest.main()
