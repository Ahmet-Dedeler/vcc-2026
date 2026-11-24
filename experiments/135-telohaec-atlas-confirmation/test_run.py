#!/usr/bin/env python3
"""Synthetic tests for Experiment 135. No TeloHAEC expression is opened."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp135_run", HERE / "run.py")
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUN)


class Experiment135Tests(unittest.TestCase):
    def test_frozen_candidate_provenance_axis_and_selection(self) -> None:
        protocol = RUN.read_json(RUN.PROTOCOL)
        panel = RUN.read_axis(RUN.PANEL, "target_gene")
        residual, confidence, support, provenance = RUN.load_candidate(
            panel, protocol
        )
        self.assertEqual(residual.shape, (300, 300))
        self.assertEqual(confidence.shape, (300, 300))
        self.assertEqual(support.shape, (300, 300))
        self.assertTrue(all(provenance.values()))

    def test_response_bins_are_balanced_and_stable(self) -> None:
        values = np.arange(300, dtype=np.float64)
        bins = RUN.response_bins(values)
        self.assertEqual(set(map(int, bins)), set(range(10)))
        self.assertTrue(
            all(np.count_nonzero(bins == index) == 30 for index in range(10))
        )

    def test_null_banks_preserve_structural_self_zero(self) -> None:
        bins = RUN.response_bins(np.arange(300, dtype=np.float64))
        target_indices = np.asarray([2, 11, 47, 101, 222], dtype=np.int16)
        target, target_response, response = RUN.freeze_nulls(
            target_indices,
            bins,
            7,
            target_seed=1352026,
            response_seed=1352027,
        )
        self.assertEqual(target.shape, (7, 5))
        self.assertEqual(target_response.shape, (7, 300))
        self.assertEqual(response.shape, (7, 5, 300))
        matrix = np.arange(300 * 300, dtype=np.float64).reshape(300, 300)
        np.fill_diagonal(matrix, 0.0)
        for draw in range(7):
            self.assertEqual(set(map(int, target[draw])), set(range(5)))
            for row, destination in enumerate(target_indices):
                source = target_indices[target[draw, row]]
                moved = matrix[source][target_response[draw]]
                self.assertEqual(moved[destination], 0.0)
                row_mapping = response[draw, row]
                self.assertEqual(row_mapping[destination], destination)
                for bin_id in range(10):
                    members = np.flatnonzero(bins == bin_id)
                    movable = members[members != destination]
                    self.assertEqual(
                        set(map(int, row_mapping[movable])),
                        set(map(int, movable)),
                    )

    def test_control_and_target_scans_read_disjoint_rows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tiny.h5"
            with h5py.File(path, "w") as handle:
                obs = handle.create_group("obs")
                obs.create_dataset(
                    "barcodes",
                    data=np.asarray(
                        ["x:CTRL:y", "x:TARG:y", "x:OTHER:y"],
                        dtype=h5py.string_dtype("utf-8"),
                    ),
                )
                x = handle.create_group("X")
                x.create_dataset(
                    "data",
                    data=np.asarray([2, 3, 5, 7, 11], dtype=np.float32),
                )
                x.create_dataset(
                    "indices",
                    data=np.asarray([0, 2, 0, 1, 2], dtype=np.int32),
                )
                x.create_dataset(
                    "indptr",
                    data=np.asarray([0, 2, 4, 5], dtype=np.int32),
                )
            with h5py.File(path, "r") as handle:
                sums, totals, libraries = RUN.scan_controls(
                    handle,
                    ("CTRL",),
                    {"CTRL": 1},
                    np.asarray([0, 1, 2], dtype=np.int32),
                    3,
                )
                panel, target_totals, spacers, rows = RUN.scan_targets(
                    handle,
                    ("GENE",),
                    {"GENE": ["TARG"]},
                    {"TARG": 1},
                    np.asarray([0, 1, -1], dtype=np.int16),
                )
        self.assertTrue(np.array_equal(sums[0], np.asarray([2.0, 0.0, 3.0])))
        self.assertTrue(np.array_equal(totals, np.asarray([5.0])))
        self.assertTrue(np.array_equal(libraries, np.asarray([5.0])))
        self.assertEqual(panel[0, 0], 5.0)
        self.assertEqual(panel[0, 1], 7.0)
        self.assertEqual(target_totals[0], 12.0)
        self.assertEqual(spacers, ("TARG",))
        self.assertEqual(rows, {"GENE": [0]})

    def test_target_metrics_rank_by_prediction_magnitude(self) -> None:
        truth = np.ones(300, dtype=np.float64)
        prediction = np.ones(300, dtype=np.float64) * 0.5
        credit, precision = RUN.target_metrics(prediction, truth, 17)
        self.assertAlmostEqual(credit, 0.5)
        self.assertEqual(precision, 1.0)
        prediction[:60] = -3.0
        credit, precision = RUN.target_metrics(prediction, truth, 299)
        self.assertLess(credit, 0.0)
        self.assertEqual(precision, 0.0)

    def test_panel_transport_is_exact_experiment_128_semantics(self) -> None:
        exp128 = RUN.load_exp128()
        rng = np.random.default_rng(7)
        genes = 500
        panel_columns = np.arange(300, dtype=np.int64)
        panel_mask = np.arange(genes) < 300
        incumbent_profile = rng.random(genes)
        incumbent_profile /= incumbent_profile.sum()
        b_ctrl = np.log1p(RUN.TS * incumbent_profile)
        incumbent_totals = np.floor(
            RUN.CELLS * 20000 * incumbent_profile
        ).astype(np.int64)
        residual = rng.normal(0, 0.05, 300)
        target = 23
        predicted, audit = RUN.transport_prediction(
            exp128,
            incumbent_profile,
            incumbent_totals,
            b_ctrl,
            panel_columns,
            panel_mask,
            residual,
            target,
        )
        full_candidate, full_audit = exp128.integer_budget_transport(
            incumbent_totals,
            incumbent_profile,
            b_ctrl,
            residual,
            panel_columns,
            panel_mask,
            target,
            1.0,
            ts=RUN.TS,
            min_cpm=RUN.MIN_CPM,
        )
        total = int(full_candidate.sum())
        expected = (
            np.log1p(RUN.TS * full_candidate[panel_columns] / total)
            - b_ctrl[panel_columns]
        )
        self.assertTrue(np.array_equal(predicted, expected))
        self.assertTrue(audit["all_invariants_pass"])
        self.assertTrue(full_audit["all_invariants_pass"])

    def test_finite_coordinate_is_mutable_even_without_source_support(self) -> None:
        exp128 = RUN.load_exp128()
        genes = 400
        panel_columns = np.arange(300, dtype=np.int64)
        panel_mask = np.arange(genes) < 300
        incumbent_profile = np.ones(genes, dtype=np.float64) / genes
        b_ctrl = np.log1p(RUN.TS * incumbent_profile)
        incumbent_totals = np.full(genes, 20000, dtype=np.int64)
        residual = np.zeros(300, dtype=np.float64)
        residual[11] = 0.5
        predicted, audit = RUN.transport_prediction(
            exp128,
            incumbent_profile,
            incumbent_totals,
            b_ctrl,
            panel_columns,
            panel_mask,
            residual,
            23,
        )
        self.assertNotEqual(predicted[11], 0.0)
        self.assertEqual(audit["mutable_coordinates"], 299)

    def test_replayed_libraries_are_deterministic_and_no_replacement(self) -> None:
        libraries = np.arange(1, 1001, dtype=np.float64)
        first = RUN.replay_libraries(
            libraries, 300, cells=400, seed=0, maximum=1_000_000
        )
        second = RUN.replay_libraries(
            libraries, 300, cells=400, seed=0, maximum=1_000_000
        )
        self.assertTrue(np.array_equal(first, second))
        self.assertEqual(first.shape, (300,))

    def test_incumbent_rows_are_target_specific(self) -> None:
        exp128 = RUN.load_exp128()
        effects = np.zeros((300, 500), dtype=np.float64)
        effects[2, 350] = 0.5
        effects[11, 351] = -0.5
        target_indices = np.asarray([2, 11], dtype=np.int16)
        base = np.ones(500, dtype=np.float64) / 500
        b_ctrl = np.log1p(RUN.TS * base)
        panel_columns = np.arange(300, dtype=np.int64)
        libraries = np.full(300, 20_000, dtype=np.int64)
        profiles, totals, predictions = RUN.build_incumbent_rows(
            exp128,
            effects,
            target_indices,
            b_ctrl,
            panel_columns,
            libraries,
        )
        self.assertEqual(profiles.shape, (2, 500))
        self.assertEqual(totals.shape, (2, 500))
        self.assertEqual(predictions.shape, (2, 300))
        self.assertFalse(np.array_equal(profiles[0], profiles[1]))

    def test_guide_diagnostics_include_cells_signs_and_leave_one_out(self) -> None:
        target_order = ("A",)
        target_indices = np.asarray([0], dtype=np.int16)
        target_rows = {"A": [0, 1, 2]}
        control = np.zeros(300, dtype=np.float64)
        guide_profiles = np.vstack(
            [
                np.linspace(-0.2, 0.2, 300),
                np.linspace(-0.19, 0.21, 300),
                np.linspace(-0.21, 0.19, 300),
            ]
        )
        truth, diagnostics = RUN.guide_truth_and_diagnostics(
            target_order,
            target_indices,
            target_rows,
            guide_profiles,
            control,
            {"G1": 10, "G2": 11, "G3": 12},
            ("G1", "G2", "G3"),
        )
        self.assertEqual(truth.shape, (1, 300))
        self.assertEqual(int(diagnostics["guide_counts"][0]), 3)
        self.assertEqual(int(diagnostics["target_cell_counts"][0]), 33)
        self.assertGreater(float(diagnostics["guide_sign_agreement"][0]), 0.95)
        self.assertGreater(
            float(diagnostics["held_out_guide_sign_agreement_min"][0]), 0.95
        )

    def test_add_one_p_value(self) -> None:
        null = np.asarray([-1.0, 0.0, 1.0])
        self.assertEqual(RUN.add_one_p(null, 2.0), 0.25)
        self.assertEqual(RUN.add_one_p(null, 0.0), 0.75)


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Experiment135Tests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)
