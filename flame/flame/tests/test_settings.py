""" Tests for the `settings` module. Variable handling and helper functions.
"""
from __future__ import print_function, division, generators
import unittest
import mock
import sys
from os import stat, environ, listdir
from os.path import join, dirname, abspath

from flame import settings as S


MOCK_DIR = join(dirname(abspath(__file__)), 'mock_dir')


def grab_mock():
    return listdir(MOCK_DIR)


class TestSettings(unittest.TestCase):
    def setUp(self):
        self._environ = environ.copy()
        self.testpath = '/some/path/we/choose'

    def test_output_dir_with_env(self):
        environ['FLAME_OUTPUT'] = self.testpath
        self.assertEqual(S.get_output_folder(), self.testpath)

    def test_output_dir_without_env(self):
        """ Here we first delete the environment variable to test what happens if none is
        provided. We climb from the this testfile to the desired output directory. Since
        there might be ambiguities in the pathname due to links, we check the indode of
        the folders.
        """
        del(environ['FLAME_OUTPUT'])
        fpath = join(S.parent(__file__), '../../../output')

        ref_inode = stat(fpath).st_ino
        test_inode = stat(S.get_output_folder()).st_ino
        self.assertEqual(ref_inode, test_inode)

    def test_failed_arrow_import(self):
        sys.modules['arrow'] = None
        notime = S.get_time()
        self.assertEqual(notime, 'n0timet0day')

    def test_arrow_import(self):
        sys.modules['arrow'] = mock.Mock()
        S.get_time()

    @mock.patch("builtins.open", mock.mock_open())
    def test_get_params(self):
        mock_params = S.get_params()
        self.assertIs(mock_params, None)

    def test_parameter_skeleton(self):
        test_skel = S.get_skel('test_skel')
        self.assertIsInstance(test_skel, str)

    def test_colors(self):
        """ This is just a quick check with a length one list, can be many more.
        """
        test_twins = (1,)
        ref_colors = [(0, 1, (68.08602, 1.24287, 84.000825))]
        test_cols = list(S.get_colors(test_twins))
        self.assertEqual(ref_colors, test_cols)

    def test_atomsIO(self):
        atom_test = S.AtomsIO('Carbon', (3, 5, 9))
        self.assertEqual(atom_test.element, 'Carbon')
        self.assertEqual(atom_test.location[0], 3)
        self.assertEqual(atom_test.location[1], 5)
        self.assertEqual(atom_test.location[2], 9)

    def tearDown(self):
        environ.clear()
        environ.update(self._environ)
