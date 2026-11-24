import math
import unittest

import numpy as np

from estimator329 import M, extractor_delta, geometry, replace_deltas, validate_fit


def fixture():
    y = np.arange(36, dtype=np.float32).reshape(3, 4, 3)/10 - .5
    se = np.full_like(y, .4)
    return y, se, extractor_delta(y), np.ones((3, 3), bool), np.ones(4, bool)


def fit(kind="iso"):
    specifications = M.components(diagonal=True)
    weights = np.zeros(len(specifications))
    weights[0 if kind == "zero" else 4] = 1  # tau=.4, variance=.16
    return {"converged": True, "components": specifications, "weights": weights.tolist()}


class InferenceTests(unittest.TestCase):
    def test_extractor_conversion(self):
        x = np.array([.1, 1., -5., np.nan, np.inf, -np.inf], np.float32)
        expected = np.where(np.isfinite(x), x*math.log(2.), 0).astype(np.float32)
        np.testing.assert_array_equal(extractor_delta(x), expected)
        self.assertEqual(extractor_delta(x).dtype, np.float32)

    def test_python_scalar_rounding_matches_extractor(self):
        x = np.arange(100, dtype=np.float32)/13
        np.testing.assert_array_equal(extractor_delta(x), x*math.log(2.))

    def test_isotropic_closed_form(self):
        y, s, d, present, genes = fixture()
        got, eligible = replace_deltas(d, y, s, present, genes, fit(), diagonal=True)
        expected = extractor_delta(y.astype(float)*(.16/(s.astype(float)**2+.16)))
        np.testing.assert_array_equal(got, expected)
        self.assertTrue(eligible.all())

    def test_zero_component(self):
        y, s, d, present, genes = fixture()
        got, _ = replace_deltas(d, y, s, present, genes, fit("zero"), diagonal=True)
        np.testing.assert_array_equal(got, 0)

    def test_batch_invariance(self):
        y, s, d, present, genes = fixture()
        a, _ = replace_deltas(d, y, s, present, genes, fit(), diagonal=True, batch_size=1)
        b, _ = replace_deltas(d, y, s, present, genes, fit(), diagonal=True, batch_size=100)
        np.testing.assert_array_equal(a, b)

    def test_missing_invalid_and_excluded_fallback(self):
        y, s, d, present, genes = fixture()
        present[0, 1] = False
        s[1, 1, 1] = 0
        s[1, 2, 2] = np.inf
        y[2, 2, 0] = np.nan
        genes[3] = False
        got, eligible = replace_deltas(d, y, s, present, genes, fit(), diagonal=True)
        self.assertFalse(eligible[0].any())
        self.assertFalse(eligible[:, 3].any())
        self.assertFalse(eligible[1, 1] or eligible[1, 2] or eligible[2, 2])
        np.testing.assert_array_equal(got[~eligible], d[~eligible])

    def test_shape_and_dtype_rejected(self):
        y, s, d, present, genes = fixture()
        with self.assertRaises(ValueError):
            replace_deltas(d, y[:-1], s, present, genes, fit(), diagonal=True)
        with self.assertRaises(ValueError):
            replace_deltas(d.astype(float), y, s, present, genes, fit(), diagonal=True)
        with self.assertRaises(ValueError):
            replace_deltas(d, y, s, present[:, :2], genes, fit(), diagonal=True)

    def test_fit_rejected(self):
        f = fit()
        f["converged"] = False
        with self.assertRaises(ValueError):
            validate_fit(f, diagonal=True)
        f = fit()
        f["weights"][0] = np.nan
        with self.assertRaises(ValueError):
            validate_fit(f, diagonal=True)
        with self.assertRaises(ValueError):
            validate_fit(fit(), diagonal=False)

    def test_geometry_identity_and_sign_change(self):
        a = np.array([[1., -1.], [0., 0.]])
        self.assertEqual(geometry(a, a)["relative_change_frobenius"], 0.)
        self.assertEqual(geometry(-a, a)["sign_changes_nonzero_baseline"], 2)
        self.assertEqual(geometry(a, a)["row_norm_ratio_min_median_max"], [1., 1., 1.])

    def test_derangement_reproducible(self):
        import run329 as A
        p = A.R.derangement(85, 20260909)
        np.testing.assert_array_equal(p, A.R.derangement(85, 20260909))
        self.assertTrue(np.all(p != np.arange(85)))
        self.assertEqual(sorted(p), list(range(85)))


if __name__ == "__main__":
    unittest.main()
