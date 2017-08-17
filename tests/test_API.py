import unittest

class testAPI(unittest.TestCase):
    def test_quick_sim(self):
        from fwdpy11_arg_example.evolve_arg import evolve_track_wrapper
        x = evolve_track_wrapper()

if __name__ == "__main__":
    unittest.main()
