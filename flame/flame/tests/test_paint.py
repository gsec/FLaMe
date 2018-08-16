import unittest
from numpy import float64
from os import rmdir
from os.path import isdir, join

# from flame.settings import HDF_EXT
from flame.tests.test_settings import grab_mock, MOCK_DIR
from flame import paint as P


def grab_first_mock():
    first_grab = [x for x in grab_mock() if x.startswith('mock_dir_0x')][0]
    return join(MOCK_DIR, first_grab)


class TestInfoExtractor(unittest.TestCase):
    def test_print_extracted_info(self):
        """ Grabs first file that matches the extension and outputs  to the logger.
        """
        P.extractor(grab_first_mock())


class TestColsOutput(unittest.TestCase):
    """ Test for correct output path creation.
    """
    def setUp(self):
        self.newdir, _, _ = P.cols_output('graph_mockdir')

    def test_dir_exists(self):
        self.assertTrue(isdir(self.newdir))

    def tearDown(self):
        rmdir(self.newdir)


class TestAveragedPlanes(unittest.TestCase):
    """ Arguably a too large of a function, nevertheless we do some basic checks.
    """
    def test_plane_attribs(self):
        test_gen = P.averaged_planes(grab_first_mock())
        next(test_gen)
        plane, color, meanflake = next(test_gen)

        print("\nPlane: {pl}\nColor: {cl}".format(pl=plane, cl=color))
        self.assertEqual(plane, '{0, 2}')

    def test_choices(self):
        choice = [1]
        test_gen = P.averaged_planes(grab_first_mock(), choice)

        for plane, color, meanflake in test_gen:
            print("\nPlane: {pl}\nColor: {cl}".format(pl=plane, cl=color))

        def dare_truth(nn, aa):
            for no, al in zip(nn, aa):
                if type(no) in (str, int, float, float64):
                    result = no == al
                elif type(no) is tuple:
                    result = dare_truth(no, al)
                else:
                    result = True
                self.assertTrue(result)
            return True

        allp = [0, 1]
        nochoic = P.averaged_planes(grab_first_mock())
        alchoic = P.averaged_planes(grab_first_mock(), allp)
        dare_truth(nochoic, alchoic)
