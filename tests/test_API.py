import unittest

class tests_ArgSimplifier(unittest.TestCase):
    def test_callable(self):
        from fwdpy11_arg_example.argsimplifier import ArgSimplifier
        self.assertEqual(callable(ArgSimplifier),True)
        a = ArgSimplifier(10)
        self.assertEqual(callable(a),True)

class testAPI(unittest.TestCase):
    def test_quick_sim(self):
        from fwdpy11_arg_example.evolve_arg import evolve_track_wrapper
        x = evolve_track_wrapper()
        print(len(x[0].gametes),len(x[0].mutations))

if __name__ == "__main__":
    unittest.main()
