from __future__ import annotations

import importlib.util
import math
import sys
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "experiment_118_pan_confirmation",
    HERE / "evaluate_pan_confirmation_v1.py",
)
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


class PanConfirmationTests(unittest.TestCase):
    def test_intersection_masks_self_and_repairs_exact_zero_sum(self) -> None:
        selected = np.array(
            [
                [0.0, 2.0, -1.0, 5.0],
                [3.0, 0.0, -2.0, 1.0],
                [2.0, -4.0, 0.0, 6.0],
                [7.0, 2.0, -3.0, 0.0],
            ]
        )
        measured = np.array([True, True, False, True])
        residual, mutable = RUN.intersected_panel_residual(
            selected, 0, measured
        )
        self.assertEqual(residual[0], 0.0)
        self.assertEqual(residual[2], 0.0)
        self.assertTrue(np.all(residual[~mutable] == 0))
        self.assertEqual(math.fsum(map(float, residual)), 0.0)

    def test_additive_application_preserves_pds_visible_profile_exactly(self) -> None:
        control = np.array([100.0, 25.0, 80.0, 40.0, 60.0, 120.0])
        panel_indices = np.array([1, 3, 4])
        residual = np.array([2.0, -2.0, 0.0])
        candidate, applied, invariants = RUN.apply_additive_linear_mean(
            control,
            panel_indices,
            residual,
            target_index=2,
            mutable_panel=np.array([True, True, False]),
        )
        self.assertTrue(invariants["all_pass"])
        self.assertEqual(math.fsum(map(float, applied)), 0.0)
        self.assertTrue(np.array_equal(candidate[[0, 2, 5]], control[[0, 2, 5]]))
        self.assertEqual(candidate[panel_indices[2]], control[panel_indices[2]])
        self.assertEqual(
            math.fsum(map(float, candidate[panel_indices])),
            math.fsum(map(float, control[panel_indices])),
        )

    def test_application_never_exponentiates_candidate_residual(self) -> None:
        control = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        panel_indices = np.array([1, 2, 3])
        residual = np.array([1.5, -1.5, 0.0])
        candidate, _, _ = RUN.apply_additive_linear_mean(
            control,
            panel_indices,
            residual,
            target_index=2,
            mutable_panel=np.array([True, True, False]),
        )
        self.assertEqual(candidate[1], 21.5)
        self.assertEqual(candidate[2], 28.5)
        self.assertEqual(candidate[3], 40.0)

    def test_factored_real_invariant_loop_retains_mutable_mask(self) -> None:
        selected = np.array(
            [
                [0.0, 2.0, -2.0],
                [1.0, 0.0, -1.0],
                [3.0, -3.0, 0.0],
            ]
        )
        control = np.array([70.0, 20.0, 80.0, 30.0, 40.0])
        panel_indices = np.array([1, 3, 4])
        measured = np.array(
            [
                [True, True, True],
                [True, True, True],
            ]
        )
        summary = RUN.application_invariant_summary(
            selected,
            control,
            panel_indices,
            measured,
            np.array([0, 1]),
        )
        self.assertTrue(summary["all_pass"])

    def test_target_equal_metrics_do_not_pool_large_targets(self) -> None:
        residual = np.array([2.0, -1.0, 0.5, 0.0])
        truth = np.array([20.0, -10.0, -5.0, 1.0])
        eligible = np.ones(4, dtype=bool)
        credit, precision, count = RUN.target_metrics(
            residual, truth, eligible
        )
        self.assertEqual(count, 4)
        self.assertAlmostEqual(credit, 2.5 / 36.0)
        self.assertEqual(precision, 1.0)

    def test_target_and_response_nulls_reapply_destination_self_exclusion(self) -> None:
        selected = np.array(
            [
                [0.0, 1.0, -1.0],
                [2.0, 0.0, -2.0],
                [4.0, -4.0, 0.0],
            ]
        )
        measured = np.ones(3, dtype=bool)
        target_null, _ = RUN.intersected_panel_residual(
            selected, 0, measured, source_target=1
        )
        response_null, _ = RUN.intersected_panel_residual(
            selected,
            0,
            measured,
            response_mapping=np.array([1, 2, 0]),
        )
        self.assertEqual(target_null[0], 0.0)
        self.assertEqual(response_null[0], 0.0)
        self.assertEqual(math.fsum(map(float, target_null)), 0.0)
        self.assertEqual(math.fsum(map(float, response_null)), 0.0)

    def test_promotion_requires_every_frozen_gate(self) -> None:
        nulls = {
            "target_credit": np.zeros(100),
            "target_precision": np.full(100, 0.5),
            "response_credit": np.zeros(100),
            "response_precision": np.full(100, 0.5),
        }
        report = RUN.promotion_report(
            0.01,
            0.60,
            {"lower_95": 0.001, "upper_95": 0.02},
            {"lower_95": 0.51, "upper_95": 0.70},
            nulls,
            64,
            True,
        )
        self.assertEqual(report["decision"], "GO")
        failed = RUN.promotion_report(
            0.01,
            0.60,
            {"lower_95": 0.0, "upper_95": 0.02},
            {"lower_95": 0.51, "upper_95": 0.70},
            nulls,
            64,
            True,
        )
        self.assertEqual(failed["decision"], "NO_GO")

    def test_authentication_does_not_numpy_load_pan_values(self) -> None:
        with mock.patch.object(
            RUN.np, "load", side_effect=AssertionError("np.load forbidden")
        ):
            authenticated = RUN.authenticate_pre_pan()
        self.assertEqual(
            len(authenticated["exact_target_intersection"]), RUN.EXPECTED_TARGETS
        )
        self.assertFalse(RUN.ACCESS_LOCK_PATH.exists())
        self.assertFalse(RUN.RESULT_DIRECTORY.exists())


if __name__ == "__main__":
    unittest.main()
