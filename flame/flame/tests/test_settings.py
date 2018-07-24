""" Tests for the `settings` module. Variable handling and helper functions.
"""
from __future__ import print_function, division, generators
import unittest
from mock import Mock
from sys import modules
from os import stat, path, environ

from flame import settings as S


class TestFolders(unittest.TestCase):
    def setUp(self):
        self._environ = environ.copy()
        self.testpath = '/some/path/we/choose'

    def test_output_dir_with_env(self):
        environ['FLAME_OUTPUT'] = self.testpath
        self.assertEqual(S.get_output_folder(), self.testpath)
        environ.clear()
        environ.update(self._environ)

    def test_output_dir_without_env(self):
        """ Here we first delete the environment variable to test what happens if none is
        provided. We climb from the this testfile to the desired output directory. Since
        there might be ambiguities in the pathname due to links, we check the indode of
        the folders.
        """
        del(environ['FLAME_OUTPUT'])
        fpath = path.join(S.parent(__file__), '../../../output')

        ref_inode = stat(fpath).st_ino
        test_inode = stat(S.get_output_folder()).st_ino
        self.assertEqual(ref_inode, test_inode)

    def test_failed_arrow_import(self):
        modules['arrow'] = None
        notime = S.get_time()
        self.assertEqual(notime, 'n0timet0day')

    def test_arrow_import(self):
        modules['arrow'] = Mock()
        S.get_time()


    def tearDown(self):
        environ.clear()
        environ.update(self._environ)


if __name__ == '__main__':
    unittest.main()
