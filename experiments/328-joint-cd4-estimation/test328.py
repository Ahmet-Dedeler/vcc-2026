import unittest

import numpy as np
from scipy.stats import multivariate_normal

import mixture328 as M


class MixtureTests(unittest.TestCase):
    def test_component_count(self):
        self.assertEqual(len(M.components()), 99)
        self.assertEqual(len(M.components(diagonal=True)), 8)
        directions = {tuple(c["u"]) for c in M.components() if c["kind"] == "rank1"}
        self.assertEqual(len(directions), 13)

    def test_zero_and_iso(self):
        y, se = np.array([[1., -2., 3.]]), np.array([[.4, .7, 1.1]])
        _, z = M.component(y, se, {"kind": "zero", "tau": 0.})
        np.testing.assert_array_equal(z, np.zeros_like(y))
        lp, p = M.component(y, se, {"kind": "iso", "tau": .8})
        v = se[0]**2+.8**2
        self.assertAlmostEqual(lp[0], multivariate_normal.logpdf(y[0], cov=np.diag(v)), places=12)
        np.testing.assert_allclose(p, y*.8**2/(se**2+.8**2), atol=1e-14)

    def test_rank_one_against_direct_algebra(self):
        rng = np.random.default_rng(328)
        y, se = rng.normal(size=(20, 3)), rng.uniform(.1, 2, (20, 3))
        for u in ([1., 0, 0], [1/np.sqrt(3)]*3, [0, 1/np.sqrt(2), -1/np.sqrt(2)]):
            spec = {"kind": "rank1", "tau": 1.6, "u": u}
            lp, mean = M.component(y, se, spec)
            prior = 1.6**2*np.outer(u, u)
            for i in range(len(y)):
                total = np.diag(se[i]**2)+prior
                self.assertAlmostEqual(lp[i], multivariate_normal.logpdf(y[i], cov=total), places=9)
                np.testing.assert_allclose(mean[i], prior @ np.linalg.solve(total, y[i]), atol=1e-12)
            reverse = {**spec, "u": list(-np.asarray(u))}
            np.testing.assert_array_equal(M.component(y, se, reverse)[0], lp)

    def test_large_parallel_signal_is_finite(self):
        y = np.full((2, 3), 1e5)
        se = np.full_like(y, .001)
        lp, mean = M.component(y, se, {"kind": "rank1", "tau": 3.2,
                                     "u": [1/np.sqrt(3)]*3})
        self.assertTrue(np.isfinite(lp).all())
        self.assertTrue(np.isfinite(mean).all())

    def test_invalid_measurement_fails(self):
        with self.assertRaises(ValueError):
            M.component(np.ones((2, 3)), np.zeros((2, 3)), M.components()[0])

    def test_em_monotonic_simplex_and_signal(self):
        rng = np.random.default_rng(20260909)
        truth = rng.normal(0, 2, (500, 1))*np.ones((1, 3))
        y = truth+rng.normal(size=(500, 3))
        se = np.ones_like(y)
        spec = [M.components()[0], {"kind": "rank1", "tau": 3.2,
                                  "u": [1/np.sqrt(3)]*3}]
        ll = M.log_densities(y[:400], se[:400], spec)
        fit = M.fit(ll)
        self.assertTrue(fit["converged"])
        self.assertAlmostEqual(sum(fit["weights"]), 1.)
        self.assertTrue(np.all(np.diff(fit["mean_loglik_trace"]) >= -1e-9))
        held = M.log_densities(y[400:], se[400:], spec)
        self.assertGreater(M.predictive(held, fit["weights"]).mean(), held[:, 0].mean())
        estimate = M.posterior(y[400:], se[400:], spec, fit["weights"])
        self.assertLess(np.mean((estimate-truth[400:])**2), np.mean((y[400:]-truth[400:])**2))

    def test_shuffle_preserves_pairs_and_complete_support(self):
        rng = np.random.default_rng(328)
        y, se = rng.normal(size=(20, 4, 3)), rng.uniform(.1, 2, (20, 4, 3))
        y[0, 0, 1] = np.nan
        a, b = M.shuffled_targets(y, se)
        a2, b2 = M.shuffled_targets(y, se)
        np.testing.assert_array_equal(a, a2)
        np.testing.assert_array_equal(b, b2)
        np.testing.assert_array_equal(a[:, :, 0], y[:, :, 0])
        valid = np.isfinite(y).all(2)
        np.testing.assert_array_equal(np.isfinite(a).all(2), valid)
        for g in range(4):
            for c in range(3):
                self.assertEqual(sorted(zip(y[valid[:, g], g, c], se[valid[:, g], g, c])),
                                 sorted(zip(a[valid[:, g], g, c], b[valid[:, g], g, c])))


if __name__ == "__main__":
    unittest.main()
