import tempfile
from pathlib import Path
import unittest

import numpy as np

import run326 as A
from report326 import CONTRASTS, linear_contrast


class CapDoseTests(unittest.TestCase):
    def test_cap30_changes_only_exponent(self):
        self.assertEqual(A.ARMS["cap15"][1:], A.ARMS["cap30"][1:])
        self.assertEqual(A.ARMS["cap30"][0], -3.)

    def test_arm_scope_restores_imported_module(self):
        before = A.R.ARMS.copy()
        with A.configured_arms():
            self.assertEqual(A.R.ARMS["cap30"], (-3., "cap", None))
        self.assertEqual(A.R.ARMS, before)

    def test_output_scope_does_not_mutate_base(self):
        original = A.E.R
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            with A.evaluator_output(path):
                self.assertIsNot(A.E.R, original)
                self.assertEqual(A.E.R.WORK, path)
                self.assertEqual(A.R.WORK, A.BASE / "work")
                self.assertEqual((path / "reference_de.parquet").resolve(),
                                 (A.BASE / "work/reference_de.parquet").resolve())
            self.assertIs(A.E.R, original)

    def test_interaction_is_difference_of_dose_differences(self):
        values = {"cap30": np.array([5., 7.]), "cap15": np.array([2., 3.]),
                  "plain30": np.array([4., 9.]), "plain15": np.array([3., 2.])}
        got = linear_contrast(values, CONTRASTS["cap-dose-minus-plain-dose"])
        np.testing.assert_array_equal(got, np.array([2., -3.]))

    def test_hash_mismatch_fails_closed(self):
        with self.assertRaises(AssertionError):
            A.require_same("fixture", "a", "b")


if __name__ == "__main__":
    unittest.main()
