from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("experiment127", HERE / "run.py")
assert SPEC and SPEC.loader
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


class TestExperiment127(unittest.TestCase):
    def test_nmae_perfect_panel(self) -> None:
        truth = np.array([1.0, -2.0, 3.0, -4.0])
        pred = np.zeros(4)
        sig = np.arange(4)
        self.assertIsNone(RUN.nmae(pred, truth, sig))
        sig = np.tile(np.arange(4), 3)[:10]
        truth2 = np.tile(truth, 3)
        pred2 = np.zeros_like(truth2)
        self.assertEqual(RUN.nmae(pred2, truth2, sig), 1.0)

    def test_fidelity_oracle_matches_bruteforce(self) -> None:
        value, calls = RUN.fidelity_oracle(5, 3, 8, 4)
        brute = [(3 + x) / max(5 + x, 8) for x in range(5)]
        self.assertAlmostEqual(value, max(brute))
        self.assertEqual(value, brute[calls])

    def test_deepest_reach_can_recover(self) -> None:
        matches = np.array([False] + [True] * 9)
        self.assertEqual(RUN.deepest_reach(matches), 1.0)

    def test_panel_first_reach_is_upper_bound_in_example(self) -> None:
        fixed = np.array([False, True, True, False, True])
        panel = np.ones(3, dtype=bool)
        panel_first = RUN.deepest_reach(np.concatenate([panel, fixed]))
        panel_last = RUN.deepest_reach(np.concatenate([fixed, panel]))
        self.assertGreaterEqual(panel_first, panel_last)

    def test_prediction_order_significance_first(self) -> None:
        sig = np.array([False, True, True])
        padj = np.array([0.9, 0.02, 0.01])
        p = np.array([0.8, 0.02, 0.01])
        lfc = np.array([9.0, 1.0, 1.0])
        genes = np.array(["A", "B", "C"])
        self.assertEqual(RUN.prediction_order(sig, padj, p, lfc, genes).tolist(), [2, 1, 0])


if __name__ == "__main__":
    unittest.main()
