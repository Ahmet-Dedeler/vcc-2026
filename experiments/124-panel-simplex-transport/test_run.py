#!/usr/bin/env python3
"""Synthetic tests for Experiment 124."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np


RUN_PATH = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("experiment_124_run", RUN_PATH)
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUN)


class PanelSimplexTests(unittest.TestCase):
    def test_nullable_string_h5ad_index_decoder(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "index.h5"
            with h5py.File(path, "w") as handle:
                group = handle.create_group("index")
                group.attrs["encoding-type"] = "nullable-string-array"
                group.create_dataset("mask", data=np.asarray([False, False]))
                string_type = h5py.string_dtype(encoding="utf-8")
                group.create_dataset(
                    "values",
                    data=np.asarray(["A", "B"], dtype=object),
                    dtype=string_type,
                )
            with h5py.File(path, "r") as handle:
                self.assertEqual(RUN.read_h5ad_index(handle["index"]), ("A", "B"))

    def test_naive_panel_effect_changes_nonpanel_after_normalization(self) -> None:
        b_ctrl = np.log1p(RUN.TS * np.asarray([0.2, 0.3, 0.5]))
        baseline = RUN.profile_from_effect(np.zeros(3), b_ctrl, min_cpm=0.0)
        candidate = RUN.profile_from_effect(
            np.asarray([0.5, 0.0, 0.0]),
            b_ctrl,
            min_cpm=0.0,
        )
        self.assertNotEqual(candidate[2], baseline[2])

    def test_transport_preserves_nonpanel_self_mass_and_positivity(self) -> None:
        baseline = np.asarray([0.10, 0.20, 0.30, 0.40], dtype=np.float64)
        b_ctrl = np.log1p(RUN.TS * baseline)
        residual = np.asarray([0.0, 0.8, -0.4], dtype=np.float64)
        candidate, audit = RUN.panel_simplex_transport(
            baseline,
            b_ctrl,
            residual,
            np.asarray([0, 1, 2]),
            0,
            1.0,
            min_cpm=0.0,
        )
        self.assertTrue(all(audit["invariants"].values()))
        self.assertEqual(candidate[3], baseline[3])
        self.assertEqual(candidate[0], baseline[0])
        self.assertTrue(np.all(candidate >= 0.0))
        self.assertEqual(
            RUN.math.fsum(map(float, candidate[:3])),
            RUN.math.fsum(map(float, baseline[:3])),
        )

    def test_transport_respects_expression_gate(self) -> None:
        baseline = np.asarray([0.10, 0.0, 0.40, 0.50], dtype=np.float64)
        b_ctrl = np.log1p(RUN.TS * baseline)
        residual = np.asarray([0.0, 5.0, -0.5], dtype=np.float64)
        candidate, audit = RUN.panel_simplex_transport(
            baseline,
            b_ctrl,
            residual,
            np.asarray([0, 1, 2]),
            0,
            1.0,
            min_cpm=10.0,
        )
        self.assertEqual(candidate[1], 0.0)
        self.assertTrue(all(audit["invariants"].values()))

    def test_full_panel_compensated_sum_is_exact(self) -> None:
        panel_size = 300
        raw = np.linspace(1.0, 3.0, panel_size + 1, dtype=np.float64)
        baseline = raw / raw.sum()
        b_ctrl = np.log1p(RUN.TS * baseline)
        residual = np.linspace(-0.9, 0.7, panel_size, dtype=np.float64)
        residual[1] = 0.0
        candidate, audit = RUN.panel_simplex_transport(
            baseline,
            b_ctrl,
            residual,
            np.arange(panel_size, dtype=np.int64),
            1,
            0.1,
            min_cpm=0.0,
        )
        self.assertTrue(all(audit["invariants"].values()))
        self.assertEqual(
            RUN.math.fsum(map(float, candidate[:panel_size])),
            RUN.math.fsum(map(float, baseline[:panel_size])),
        )

    def test_emitted_nonpanel_totals_are_exact_when_profiles_are_exact(self) -> None:
        baseline = np.asarray([0.1, 0.2, 0.3, 0.4], dtype=np.float64)
        candidate = baseline.copy()
        candidate[0], candidate[1] = 0.15, 0.15
        nonpanel = np.asarray([False, False, True, True])
        self.assertTrue(
            RUN.emitted_nonpanel_totals_identical(
                baseline,
                candidate,
                nonpanel,
                (19_946, 20_034, 20_109),
                400,
            )
        )


if __name__ == "__main__":
    unittest.main()
