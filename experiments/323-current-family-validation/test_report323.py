import unittest
import numpy as np

from report323 import mse_bootstrap
from cell_eval2.metrics.delta import _across_pert_budget


class MseBootstrapTests(unittest.TestCase):
    def test_weighted_budget_matches_direct_duplicate_rows(self):
        rng = np.random.default_rng(323)
        p = rng.normal(size=(4, 9))
        r = rng.normal(size=(4, 9))
        own = np.arange(4)
        draws = np.array([[1, 1, 1, 1], [0, 2, 1, 1], [4, 0, 0, 0]])
        _, budget, _, _ = mse_bootstrap(
            draws, p, r, np.zeros(9), np.full(4, .3), np.full(4, .4), .1, own)
        for i, counts in enumerate(draws):
            idx = np.repeat(np.arange(4), counts)
            direct = _across_pert_budget(
                p[idx], {j: int(own[k]) for j, k in enumerate(idx)},
                np.full(4, 1/8))
            self.assertAlmostEqual(budget[i], direct, places=12)


if __name__ == "__main__":
    unittest.main()
