import unittest
import os

# from flame.settings import GRAPH_OUTPUT
from flame import paint as P


class TestColsOutput(unittest.TestCase):
    """ Test for correct output path creation.
    """
    def setUp(self):
        self.newdir, _, _ = P.cols_output('mockdirfortesting_creat')

    def test_dir_exists(self):
        self.assertTrue(os.path.isdir(self.newdir))

    def tearDown(self):
        os.rmdir(self.newdir)


if __name__ == '__main__':
    unittest.main()
