"""Small fail-closed tests, with no donor or reference data access."""
import unittest

import numpy as np

from run323 import axis_columns, axis_digest, derangement, map_axis, verify_axis


class AxisTests(unittest.TestCase):
    def test_panel_excluded(self):
        genes = np.array(["A", "B", "C", "D"])
        cols = axis_columns(genes, [["A", "B"], ["B", "C"]], ["B"])
        np.testing.assert_array_equal(cols, [0, 2])
        verify_axis(genes, cols, ["B"], axis_digest(["A", "C"]))

    def test_wrong_panel_rejected(self):
        with self.assertRaisesRegex(ValueError, "excluded"):
            verify_axis(np.array(["A", "B", "C"]), [0, 1], ["B"],
                        axis_digest(["A", "B"]))

    def test_reordered_axis_rejected(self):
        with self.assertRaisesRegex(ValueError, "digest"):
            verify_axis(np.array(["A", "B", "C"]), [1, 0], ["C"],
                        axis_digest(["A", "B"]))

    def test_mapping_by_symbol_and_missing_zero(self):
        np.testing.assert_array_equal(
            map_axis(np.array([[1, 2, 3]]), ["A", "B", "C"], ["C", "Z", "A"]),
            [[3, 0, 1]])

    def test_duplicate_axis_rejected(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            verify_axis(np.array(["A", "B"]), [0, 0], [], axis_digest(["A", "A"]))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            map_axis(np.ones((1, 2)), ["A", "A"], ["A"])

    def test_derangement_preserves_rows(self):
        p = derangement(85)
        self.assertTrue(np.all(p != np.arange(85)))
        np.testing.assert_array_equal(np.sort(p), np.arange(85))
        np.testing.assert_array_equal(p, derangement(85))

    def test_projected_transform_no_target_crosstalk(self):
        effect = np.array([[0., 7., 0.]])
        cols = axis_columns(["A", "B", "C"], [["A", "B", "C"]], ["B"])
        out = effect.copy()
        out[:, cols] = effect[:, cols] @ np.array([[2., 1.], [1., 3.]])
        np.testing.assert_array_equal(out[:, cols], 0)


if __name__ == "__main__":
    unittest.main()
