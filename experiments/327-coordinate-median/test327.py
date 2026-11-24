import unittest

import numpy as np

from estimator327 import arms, eligible_mask, interpolate, restore_subvector


class EstimatorTests(unittest.TestCase):
    def setUp(self):
        self.sources = np.array([[[1., 2., 3., 5.]],
                                 [[1., 3., 4., 6.]],
                                 [[-2.2, 4., 5., 7.]]])
        self.base = self.sources.mean(0)
        self.observed = np.ones_like(self.sources, dtype=bool)
        self.op = np.array([True, True, True, False])

    def test_zero_dose_is_exact(self):
        mask = eligible_mask(self.observed, self.op)
        np.testing.assert_array_equal(
            interpolate(self.base, np.median(self.sources, 0), mask, 0), self.base)

    def test_median_can_repair_single_outlier_sign(self):
        got, _ = arms(self.base, self.sources, self.observed, self.op)
        self.assertLess(self.base[0, 0], 0)
        self.assertGreater(got["median25"][0, 0], 0)

    def test_missing_coordinate_is_not_a_zero_observation(self):
        self.observed[2, 0, 0] = False
        got, mask = arms(self.base, self.sources, self.observed, self.op)
        self.assertFalse(mask[0, 0])
        for value in got.values():
            self.assertEqual(value[0, 0], self.base[0, 0])

    def test_missing_target_is_identity(self):
        self.observed[2] = False
        got, mask = arms(self.base, self.sources, self.observed, self.op)
        self.assertFalse(mask.any())
        for value in got.values():
            np.testing.assert_array_equal(value, self.base)

    def test_norm_and_outside_support(self):
        got, mask = arms(self.base, self.sources, self.observed, self.op)
        for value in got.values():
            np.testing.assert_array_equal(value[~mask], self.base[~mask])
            self.assertAlmostEqual(np.linalg.norm(value[mask]),
                                   np.linalg.norm(self.base[mask]), places=12)

    def test_sign_control_magnitudes_are_exact(self):
        got, _ = arms(self.base, self.sources, self.observed, self.op)
        np.testing.assert_array_equal(np.abs(got["magnitude25"]), np.abs(got["median25"]))
        np.testing.assert_array_equal(np.sign(got["magnitude25"]), np.sign(self.base))
        np.testing.assert_array_equal(np.sign(got["power25"]), np.sign(self.base))

    def test_shuffled_correction_is_deterministic(self):
        first, _ = arms(self.base, self.sources, self.observed, self.op)
        second, _ = arms(self.base, self.sources, self.observed, self.op)
        np.testing.assert_array_equal(first["shuffled25"], second["shuffled25"])
        self.assertFalse(np.array_equal(first["shuffled25"], first["median25"]))

    def test_zero_and_nonfinite_guard(self):
        mask = np.ones_like(self.base, bool)
        np.testing.assert_array_equal(restore_subvector(
            self.base, np.zeros_like(self.base), mask), self.base)
        with self.assertRaises(ValueError):
            restore_subvector(self.base, np.full_like(self.base, np.nan), mask)
        with self.assertRaises(ValueError):
            interpolate(self.base, self.base, mask, 2)


if __name__ == "__main__":
    unittest.main()
