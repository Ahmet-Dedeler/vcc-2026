from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("experiment_123_run", HERE / "run.py")
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


class SlimReproductionTests(unittest.TestCase):
    def test_candidate_grid_is_frozen_and_contains_identity(self) -> None:
        grid = RUN.candidates()
        self.assertEqual(len(grid), 73)
        self.assertEqual(grid[0].key, "identity")
        self.assertIn("slim-blend:K=10:ridge=0.1:f=0.1", {c.key for c in grid})

    def test_stable_pca_scores_are_deterministic_and_orthogonal(self) -> None:
        rng = np.random.default_rng(1)
        y = rng.normal(size=(40, 12))
        first = RUN.stable_pca_scores(y, 5)
        second = RUN.stable_pca_scores(y, 5)
        np.testing.assert_array_equal(first, second)
        offdiag = first.T @ first
        np.testing.assert_allclose(offdiag, np.diag(np.diag(offdiag)), atol=1e-9)

    def test_closed_form_matches_upstream_equation(self) -> None:
        rng = np.random.default_rng(2)
        genes, targets, K, dim = 25, 18, 5, 7
        G = rng.normal(size=(genes, K))
        P = rng.normal(size=(targets, dim))
        Y = rng.normal(size=(genes, targets))
        factors = {"G": G, "gt_y": G.T @ Y}
        observed = RUN.solve_w(factors, P, K, 0.1)
        z = np.linalg.solve(G.T @ G + 0.1 * np.eye(K), G.T @ Y)
        expected = np.linalg.solve(P.T @ P + 0.1 * np.eye(dim), (z @ P).T).T
        np.testing.assert_allclose(observed, expected, atol=1e-12)

    def test_direction_blend_identity_and_norm_preservation(self) -> None:
        base = np.asarray([[3.0, 4.0], [1.0, 0.0]])
        learned = np.asarray([[4.0, -3.0], [0.0, 2.0]])
        np.testing.assert_array_equal(RUN.direction_blend(base, learned, 0.0), base)
        mixed = RUN.direction_blend(base, learned, 0.35)
        np.testing.assert_allclose(
            np.linalg.norm(mixed, axis=1), np.linalg.norm(base, axis=1), atol=1e-12
        )

    def test_direction_blend_falls_back_for_zero_slim(self) -> None:
        base = np.asarray([[1.0, 2.0], [3.0, 4.0]])
        learned = np.asarray([[0.0, 0.0], [1.0, -1.0]])
        mixed = RUN.direction_blend(base, learned, 0.5)
        np.testing.assert_array_equal(mixed[0], base[0])

    def test_pds_identity_is_one(self) -> None:
        truth = np.eye(5, dtype=np.float64)
        np.testing.assert_allclose(RUN.pds_per_target(truth, truth), 1.0)

    def test_pds_zero_row_is_chance(self) -> None:
        truth = np.eye(4, dtype=np.float64)
        prediction = truth.copy()
        prediction[2] = 0.0
        self.assertEqual(RUN.pds_per_target(prediction, truth)[2], 0.5)

    def test_signed_credit_rewards_matching_signs(self) -> None:
        truth = np.asarray([[2.0, -1.0, 0.5]])
        good = np.asarray([[3.0, -2.0, 1.0]])
        bad = -good
        self.assertGreater(
            RUN.signed_credit_per_target(good, truth)[0],
            RUN.signed_credit_per_target(bad, truth)[0],
        )

    def test_fast_blend_pds_matches_materialized_prediction(self) -> None:
        rng = np.random.default_rng(7)
        genes, fit_count, query_count, K, dim = 31, 30, 6, 5, 8
        G = rng.normal(size=(genes, RUN.MAX_K))
        y = rng.normal(size=(genes, fit_count))
        factors = {
            "G": G,
            "gt_y": G.T @ y,
        }
        fit_embeddings = rng.normal(size=(fit_count, dim))
        query_embeddings = rng.normal(size=(query_count, dim))
        incumbent = rng.normal(size=(query_count, genes))
        truth = rng.normal(size=(query_count, genes))
        fast = RUN.candidate_pds_deltas(
            factors,
            fit_embeddings,
            query_embeddings,
            incumbent,
            truth,
        )
        W = RUN.solve_w(factors, fit_embeddings, K, 0.1)
        slim = query_embeddings @ W.T @ G[:, :K].T
        candidate = RUN.direction_blend(incumbent, slim, 0.35)
        expected = float(
            np.mean(
                RUN.pds_per_target(candidate, truth)
                - RUN.pds_per_target(incumbent, truth)
            )
        )
        self.assertAlmostEqual(
            fast["slim-blend:K=5:ridge=0.1:f=0.35"], expected, places=12
        )

    def test_top20_precision_uses_prediction_magnitude(self) -> None:
        truth = np.asarray([[1.0, -5.0, 1.0, 1.0, 1.0]])
        good = np.asarray([[0.0, -10.0, 0.0, 0.0, 0.0]])
        bad = -good
        self.assertEqual(RUN.top20_precision_per_target(good, truth)[0], 1.0)
        self.assertEqual(RUN.top20_precision_per_target(bad, truth)[0], 0.0)

    def test_selection_keeps_identity_when_every_candidate_loses(self) -> None:
        deltas = {
            "a_to_b": {
                candidate.key: (0.0 if candidate.kind == "identity" else -0.1)
                for candidate in RUN.candidates()
            },
            "b_to_a": {
                candidate.key: (0.0 if candidate.kind == "identity" else -0.2)
                for candidate in RUN.candidates()
            },
        }
        selected, _ = RUN.select_candidate_deltas(deltas)
        self.assertEqual(selected.kind, "identity")

    def test_target_order_is_deterministic(self) -> None:
        targets = {"A", "B", "C", "D"}
        self.assertEqual(RUN.stable_target_order(targets), RUN.stable_target_order(reversed(sorted(targets))))

    def test_only_three_source_donors_are_executable(self) -> None:
        self.assertEqual(
            RUN.DONORS,
            ("k562_gwps", "rpe1", "nadig_jurkat"),
        )
        self.assertEqual(
            RUN.HELD_FORBIDDEN,
            ("h1_2025", "orion_hct116", "orion_hek293t"),
        )


if __name__ == "__main__":
    unittest.main()
