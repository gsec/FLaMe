import unittest
from flame.simulation import builder


class TestBuilder(unittest.TestCase):
    def test_geometry(self):
        gen_flake = builder((-2, 2), total_size=100, snapshot_interval=50)
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


if __name__ == '__main__':
    unittest.main()
