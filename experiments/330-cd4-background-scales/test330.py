from types import SimpleNamespace
import unittest

import numpy as np

import scales330 as S


class ScaleTests(unittest.TestCase):
    def test_exact_formula(self):
        x = np.arange(40).reshape(10, 4)*np.arange(1, 5)
        a = np.abs(x-x.mean(0)).mean(0)
        expected = np.maximum(a, np.percentile(a[a > 0], 25))
        got, info = S.gene_scale(x)
        np.testing.assert_array_equal(got, expected)
        self.assertEqual(info["n"], 10)

    def test_translation_invariance(self):
        rng = np.random.default_rng(330)
        x = rng.normal(size=(20, 8))
        a, _ = S.gene_scale(x)
        b, _ = S.gene_scale(x + np.arange(8)[None, :])
        np.testing.assert_allclose(a, b, atol=1e-15, rtol=1e-15)

    def test_degenerate_rejected(self):
        with self.assertRaises(ValueError):
            S.gene_scale(np.ones((10, 4)))
        with self.assertRaises(ValueError):
            S.gene_scale(np.full((10, 4), np.nan))

    def test_partial_override(self):
        got = S.mixed_scale([1, 2, 3], [2, 6, 9], [5, 7], [9, 10])
        np.testing.assert_array_equal(got, [1, 2, 5])
        with self.assertRaises(ValueError):
            S.mixed_scale([1, 2], [2, 2], [5], [2])

    def test_median_and_chunk_invariance(self):
        x = np.arange(60).reshape(10, 6)
        cols = np.array([1, 3, 5])
        mean = x[:, cols].mean(0)
        scale = np.array([1., 2., 3.])
        expected = np.median(np.linalg.norm((x[:, cols]-mean)/scale, axis=1))
        for chunk in (1, 3, 128):
            self.assertEqual(S.median_normalized_norm(x, cols, mean, scale, chunk), expected)

    def test_shuffle_marginals_and_determinism(self):
        x = np.arange(60).reshape(3, 20)+1
        a, b = S.permute_scales(x), S.permute_scales(x)
        np.testing.assert_array_equal(a, b)
        np.testing.assert_array_equal(np.sort(a, axis=1), x)
        self.assertFalse(np.array_equal(a, x))

    def test_identity_override_restores_module(self):
        def f(*args, **kwargs):
            return 2, np.array([3., 4.]), 5
        m = SimpleNamespace(centered_statistics=f)
        with S.scale_override(m, None, np.array([]), {}):
            got = m.centered_statistics(None, None, None, need_gene_scale=True,
                                       chunk_rows=128, max_rss_gib=22)
            np.testing.assert_array_equal(got[1], [3, 4])
            self.assertEqual(got[0], 2)
            self.assertEqual(got[2], 5)
        self.assertIs(m.centered_statistics, f)
        with self.assertRaises(RuntimeError):
            with S.scale_override(m, None, np.array([]), {}):
                raise RuntimeError("synthetic")
        self.assertIs(m.centered_statistics, f)


if __name__ == "__main__":
    unittest.main()
