import unittest

class tests_ArgSimplifier(unittest.TestCase):
    def test_callable(self):
        from fwdpy11_arg_example.argsimplifier import ArgSimplifier
        self.assertEqual(callable(ArgSimplifier),True)
        a = ArgSimplifier(10)
        self.assertEqual(callable(a),True)


if __name__ == "__main__":
    unittest.main()
