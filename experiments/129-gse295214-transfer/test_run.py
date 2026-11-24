from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("experiment_129_run", HERE / "run.py")
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


class Experiment129Tests(unittest.TestCase):
    def test_candidate_grid_is_low_capacity_and_contains_identity(self) -> None:
        self.assertEqual(
            [candidate.key for candidate in RUN.candidates()],
            [
                "identity",
                "pan-knn:k=5:angle=5deg",
                "pan-knn:k=15:angle=5deg",
                "pan-knn:k=30:angle=5deg",
            ],
        )

    def test_stable_order_and_folds_are_deterministic(self) -> None:
        targets = tuple(f"T{index:03d}" for index in range(100))
        self.assertEqual(
            RUN.stable_order(targets, RUN.TARGET_NAMESPACE),
            RUN.stable_order(reversed(targets), RUN.TARGET_NAMESPACE),
        )
        np.testing.assert_array_equal(
            RUN.deterministic_folds(targets),
            RUN.deterministic_folds(targets),
        )

    def test_derangement_has_no_fixed_points(self) -> None:
        permutation = RUN.derangement(np.random.default_rng(1), 55)
        self.assertEqual(sorted(permutation.tolist()), list(range(55)))
        self.assertTrue(np.all(permutation != np.arange(55)))

    def test_duplicate_gene_axis_uses_last_occurrence(self) -> None:
        bank = RUN.Bank("synthetic", Path("."), ("T",), ("A", "B", "A"))
        rows, columns = RUN.bank_indices(bank, ("T",), ("A", "B"))
        np.testing.assert_array_equal(rows, [0])
        np.testing.assert_array_equal(columns, [2, 1])

    def test_pan_row_mask_preserves_supported_values_and_exact_zero_fallback(self) -> None:
        values = np.asarray([[1.5, 0.0, -2.0], [0.0, 3.0, 0.0]])
        support = np.asarray([[True, False, True], [False, True, False]])
        masked = RUN.apply_measured_mask(values, support)
        np.testing.assert_array_equal(masked, values)
        np.testing.assert_array_equal(masked[~support], 0.0)

    def test_pan_row_mask_rejects_nonzero_unsupported_delta(self) -> None:
        values = np.asarray([[1.5, 0.25, -2.0]])
        support = np.asarray([[True, False, True]])
        with self.assertRaisesRegex(RUN.ContractError, "unsupported delta"):
            RUN.apply_measured_mask(values, support)

    def test_tangent_residual_is_orthogonal_to_incumbent(self) -> None:
        incumbent = np.asarray([[3.0, 4.0, 0.0], [1.0, 0.0, 0.0]])
        truth = np.asarray([[1.0, -2.0, 3.0], [2.0, 1.0, 0.0]])
        residual = RUN.tangent_residual(incumbent, truth)
        dot = np.sum(RUN.row_unit(incumbent) * residual, axis=1)
        np.testing.assert_allclose(dot, 0.0, atol=1e-12)

    def test_knn_rotation_preserves_norm_and_fallback(self) -> None:
        incumbent = np.asarray([[3.0, 4.0], [2.0, 0.0]])
        query_pan = np.asarray([[1.0, 0.0], [-1.0, 0.0]])
        train_pan = np.asarray([[1.0, 0.0], [0.8, 0.2]])
        residual = np.asarray([[0.0, 1.0], [0.0, 2.0]])
        prediction, used = RUN.knn_tangent_prediction(
            incumbent, query_pan, train_pan, residual, 2
        )
        self.assertTrue(used[0])
        self.assertFalse(used[1])
        np.testing.assert_allclose(
            np.linalg.norm(prediction[0]), np.linalg.norm(incumbent[0]), atol=1e-12
        )
        np.testing.assert_array_equal(prediction[1], incumbent[1])

    def test_pds_identity_is_one_and_zero_row_is_chance(self) -> None:
        truth = np.eye(5, dtype=np.float64)
        np.testing.assert_allclose(RUN.pds_per_target(truth, truth), 1.0)
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

    def test_cross_context_teacher_isolation(self) -> None:
        rng = np.random.default_rng(7)
        n, genes = 20, 8
        baseline = rng.normal(size=(n, genes))
        pan = rng.normal(size=(n, genes))
        truths = {
            "rpe1": rng.normal(size=(n, genes)),
            "nadig_jurkat": rng.normal(size=(n, genes)),
        }
        folds = np.arange(n) % RUN.FOLDS
        candidate = RUN.Candidate("knn", 5)
        original = {
            name: values.copy()
            for name, values in truths.items()
        }
        first, _ = RUN.selection_predictions(baseline, pan, truths, folds, candidate)
        truths["rpe1"] *= -3.0
        second, _ = RUN.selection_predictions(baseline, pan, truths, folds, candidate)
        np.testing.assert_array_equal(first["rpe1"], second["rpe1"])
        self.assertFalse(np.array_equal(first["nadig_jurkat"], second["nadig_jurkat"]))
        np.testing.assert_array_equal(original["nadig_jurkat"], truths["nadig_jurkat"])

    def test_identity_wins_when_nonidentity_loses(self) -> None:
        rng = np.random.default_rng(9)
        n, genes = 100, 10
        baseline = rng.normal(size=(n, genes))
        pan = rng.normal(size=(n, genes))
        truths = {
            "rpe1": baseline.copy(),
            "nadig_jurkat": baseline.copy(),
        }
        folds = np.arange(n) % RUN.FOLDS
        selected, _ = RUN.select_candidate(baseline, pan, truths, folds)
        self.assertEqual(selected.kind, "identity")

    def test_protocol_forbids_direct_pan_and_held_access(self) -> None:
        protocol = RUN.read_json(RUN.PROTOCOL_PATH)
        self.assertIn("direct_same_target_Pan_residual", protocol["correction_boundary"]["closed_families"])
        self.assertFalse(protocol["authorization_boundary"]["held_outcome_access"])
        self.assertFalse(protocol["authorization_boundary"]["cell_eval2"])
        self.assertFalse(protocol["authorization_boundary"]["count_generation"])


if __name__ == "__main__":
    unittest.main()
