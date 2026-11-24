#!/usr/bin/env python3
"""Source-free synthetic tests for Experiment 128."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np


RUN_PATH = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("experiment_128_run", RUN_PATH)
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUN)


class IntegerBudgetTransportTests(unittest.TestCase):
    def test_source_coordinate_formula_is_exact(self) -> None:
        panel = np.asarray([0.1, 0.2, 0.3], dtype=np.float64)
        residual = np.asarray([0.0, 0.4, -0.2], dtype=np.float64)
        mutable = np.asarray([False, True, True])
        weights, negative = RUN.source_coordinate_transform(
            panel,
            residual,
            mutable,
            0.5,
            ts=50_000.0,
        )
        expected = np.expm1(
            np.log1p(50_000.0 * panel[mutable]) + 0.5 * residual[mutable]
        ) / 50_000.0
        np.testing.assert_array_equal(weights[mutable], expected)
        self.assertFalse(np.any(negative))

    def test_source_units_are_not_plain_multiplicative_fold_change(self) -> None:
        panel = np.asarray([1e-5, 0.2], dtype=np.float64)
        residual = np.asarray([0.5, 0.5], dtype=np.float64)
        mutable = np.asarray([True, True])
        weights, _ = RUN.source_coordinate_transform(
            panel,
            residual,
            mutable,
            1.0,
            ts=50_000.0,
        )
        self.assertFalse(np.array_equal(weights, panel * np.exp(residual)))

    def test_preclip_negative_is_recorded_and_clipped(self) -> None:
        panel = np.asarray([1e-5, 0.2], dtype=np.float64)
        residual = np.asarray([-2.0, 0.0], dtype=np.float64)
        mutable = np.asarray([True, False])
        weights, negative = RUN.source_coordinate_transform(
            panel,
            residual,
            mutable,
            1.0,
            ts=50_000.0,
        )
        self.assertTrue(negative[0])
        self.assertEqual(weights[0], 0.0)
        self.assertEqual(weights[1], panel[1])

    def test_largest_remainder_uses_canonical_ties(self) -> None:
        allocated, fallback = RUN.largest_remainder(
            2,
            np.ones(4, dtype=np.float64),
            np.arange(4, dtype=np.int64),
            fallback=np.asarray([1, 1, 0, 0], dtype=np.int64),
        )
        np.testing.assert_array_equal(allocated, np.asarray([1, 1, 0, 0]))
        self.assertFalse(fallback)

    def test_largest_remainder_zero_weight_fallback(self) -> None:
        expected = np.asarray([2, 0, 3], dtype=np.int64)
        allocated, fallback = RUN.largest_remainder(
            5,
            np.zeros(3, dtype=np.float64),
            np.arange(3, dtype=np.int64),
            fallback=expected,
        )
        np.testing.assert_array_equal(allocated, expected)
        self.assertTrue(fallback)

    def test_integer_transport_preserves_scorer_visible_coordinates(self) -> None:
        profile = np.asarray([0.10, 0.20, 0.30, 0.40], dtype=np.float64)
        b_ctrl = np.log1p(50_000.0 * profile)
        residual = np.asarray([0.0, 0.8, -0.4], dtype=np.float64)
        panel_columns = np.asarray([0, 1, 2], dtype=np.int64)
        panel_mask = np.asarray([True, True, True, False])
        baseline_totals = np.floor(400 * 20_000 * profile).astype(np.int64)
        candidate, audit = RUN.integer_budget_transport(
            baseline_totals,
            profile,
            b_ctrl,
            residual,
            panel_columns,
            panel_mask,
            0,
            1.0,
            ts=50_000.0,
            min_cpm=0.0,
        )
        self.assertTrue(audit["all_invariants_pass"])
        self.assertEqual(candidate[3], baseline_totals[3])
        self.assertEqual(candidate[0], baseline_totals[0])
        self.assertEqual(candidate.sum(), baseline_totals.sum())
        self.assertGreater(audit["changed_mutable_coordinates"], 0)

    def test_ineligible_panel_coordinate_is_copied(self) -> None:
        profile = np.asarray([0.10, 0.0, 0.30, 0.60], dtype=np.float64)
        b_ctrl = np.log1p(50_000.0 * profile)
        residual = np.asarray([0.0, 5.0, -0.5], dtype=np.float64)
        panel_columns = np.asarray([0, 1, 2], dtype=np.int64)
        panel_mask = np.asarray([True, True, True, False])
        baseline_totals = np.floor(400 * 20_000 * profile).astype(np.int64)
        candidate, audit = RUN.integer_budget_transport(
            baseline_totals,
            profile,
            b_ctrl,
            residual,
            panel_columns,
            panel_mask,
            0,
            1.0,
            ts=50_000.0,
            min_cpm=10.0,
        )
        self.assertTrue(audit["all_invariants_pass"])
        self.assertEqual(candidate[1], baseline_totals[1])

    def test_integer_budget_closes_flooring_denominator_leak(self) -> None:
        profile = np.asarray([0.101, 0.202, 0.303, 0.394], dtype=np.float64)
        b_ctrl = np.log1p(50_000.0 * profile)
        residual = np.asarray([0.0, 0.37, -0.19], dtype=np.float64)
        panel_columns = np.asarray([0, 1, 2], dtype=np.int64)
        panel_mask = np.asarray([True, True, True, False])
        baseline_totals = np.floor(400 * 19_946 * profile).astype(np.int64)
        candidate, audit = RUN.integer_budget_transport(
            baseline_totals,
            profile,
            b_ctrl,
            residual,
            panel_columns,
            panel_mask,
            0,
            1.0,
            ts=50_000.0,
            min_cpm=0.0,
        )
        self.assertTrue(audit["all_invariants_pass"])
        baseline_log = np.log1p(
            50_000.0 * baseline_totals / baseline_totals.sum()
        )
        candidate_log = np.log1p(candidate * 50_000.0 / candidate.sum())
        np.testing.assert_array_equal(candidate_log[3:], baseline_log[3:])


if __name__ == "__main__":
    unittest.main()
