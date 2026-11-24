import unittest
import numpy as np
from convert_canonical333 import exact_subset


class CanonicalTests(unittest.TestCase):
    def test_exact_case_whitespace_and_aliases_never_merged(self):
        source = {"C1orf1": [1., 2.], "C1ORF1": [8., 9.],
                  " ABC": [4., 5.], "xyz": [6., 7.], "CD4": [2., 3.]}
        names, matrix = exact_subset(source, ["C1orf1", "ABC", "XYZ", "CD4"])
        self.assertEqual(names.tolist(), ["C1orf1", "CD4"])
        np.testing.assert_array_equal(matrix, [[1., 2.], [2., 3.]])

    def test_unused_alias_is_not_certified(self):
        names, _ = exact_subset({"ABC": [1., 2.], "junk": [np.nan]}, ["ABC"])
        self.assertEqual(names.tolist(), ["ABC"])

    def test_selected_invalid_vector_fails(self):
        with self.assertRaises(ValueError):
            exact_subset({"ABC": [np.nan]}, ["ABC"])
        with self.assertRaises(ValueError):
            exact_subset({"ABC": [0., 0.]}, ["ABC"])

    def test_duplicate_canonical_and_dimension_mismatch_fail(self):
        with self.assertRaises(ValueError):
            exact_subset({"ABC": [1., 2.]}, ["ABC", "ABC"])
        with self.assertRaises(ValueError):
            exact_subset({"ABC": [1., 2.], "XYZ": [1.]}, ["ABC", "XYZ"])


if __name__ == "__main__":
    unittest.main()
