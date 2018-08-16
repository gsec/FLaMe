import unittest
from os import chdir
from flame import simulation as S
from flame.tests.test_settings import MOCK_DIR, grab_mock

chdir(MOCK_DIR)


class TestBuilder(unittest.TestCase):
    def test_geometry(self):
        gen_flake = S.builder((-2, 2), total_size=100, snapshot_interval=50)
        ref = (dict([
                    ('radius', 6.0), ('height', 6.531972647421807),
                    ('aspect_ratio', 1.8371173070873839), ('area', 113.09733552923255),
                    ('layers', 4), ('iter', 51)]),
               dict([
                    ('radius', 8.717797887081348), ('height', 8.246211251235321),
                    ('aspect_ratio', 2.1143765594836976), ('area', 238.76104167282432),
                    ('layers', 5), ('iter', 101)])
               )

        for testrun, reference in zip(gen_flake, ref):
            for k, v in reference.items():
                self.assertEqual(type(v), type(testrun.get(k)))


class TestSimulationRun(unittest.TestCase):
    """ Run complete Simulation from mock directory.
    """
    def setUp(self):
        self.params = {'name': 'zz_automated_testrun',
                       'function': '(0, x)',
                       'values': [1],
                       'sample_size': 2,
                       'snapshot_interval': 10,
                       'total_size': 30}

    def test_run_with_params(self):
        start_list = grab_mock()
        S.run(self.params)

        # Here we only test if the directory content changed
        self.assertNotEqual(start_list, grab_mock())
