import unittest

import numpy as np

import run331 as A


class SyntheticTests(unittest.TestCase):
    def test_shape_finite_and_reproducibility(self):
        a, b = A.synthetic(512, 331), A.synthetic(512, 331)
        self.assertEqual([x.shape for x in a], [(512, 12), (512,), (256, 12), (256,)])
        for x, y in zip(a, b):
            self.assertTrue(np.isfinite(x).all())
            np.testing.assert_array_equal(x, y)

    def test_seed_changes_problem(self):
        self.assertFalse(np.array_equal(A.synthetic(512, 331)[0], A.synthetic(512, 332)[0]))


if __name__ == "__main__":
    unittest.main()
