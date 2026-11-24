#!/usr/bin/env python3
"""Synthetic and fail-closed tests for Experiment 130."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from unittest import mock
import unittest

import numpy as np


RUN_PATH = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("experiment_130_run", RUN_PATH)
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUN)


class Experiment130Tests(unittest.TestCase):
    def test_value_access_fails_before_h5_open(self) -> None:
        self.assertFalse(RUN.BINDING.exists())
        with mock.patch.object(RUN.h5py, "File") as opener:
            with self.assertRaises(RUN.ValueAccessLocked):
                RUN.open_telohaec_values()
            opener.assert_not_called()

    def test_binding_validation_fails_before_candidate_array_load(self) -> None:
        with mock.patch.object(RUN.np, "load") as loader:
            with self.assertRaises(RUN.ValueAccessLocked):
                RUN.validate_candidate_binding()
            loader.assert_not_called()

    def test_metadata_structure_reader_never_slices_x(self) -> None:
        structure = RUN.h5ad_structure(RUN.TELO_H5AD)
        self.assertEqual(structure["x"]["shape"], [214449, 17472])
        self.assertFalse(structure["layers_present"])

    def test_guide_mapping_has_strict_replication(self) -> None:
        panel = RUN.csv_axis(RUN.PANEL, "target_gene")
        mapping = RUN.guide_mapping(panel)
        self.assertEqual(len(mapping["targets"]), 65)
        self.assertEqual(len(mapping["controls"]), 549)
        self.assertTrue(all(len(guides) >= 2 for guides in mapping["targets"].values()))

    def test_target_metrics_match_definition(self) -> None:
        truth = np.asarray([2.0, -1.0, 0.5, -0.5, 1.0])
        pred = np.asarray([1.0, -2.0, -1.0, 0.0, 0.5])
        confidence = np.asarray([5.0, 4.0, 3.0, 2.0, 1.0])
        with self.assertRaises(RUN.ExperimentError):
            RUN.target_metrics(pred, confidence, truth)
        truth = np.tile(truth, 4)
        pred = np.tile(pred, 4)
        confidence = np.arange(20, 0, -1, dtype=float)
        credit, precision = RUN.target_metrics(pred, confidence, truth)
        self.assertTrue(np.isfinite(credit))
        self.assertGreaterEqual(precision, 0.0)
        self.assertLessEqual(precision, 1.0)

    def test_derangement_has_no_fixed_points(self) -> None:
        draw = RUN.derangement(np.random.default_rng(1302026), 65)
        self.assertTrue(np.all(draw != np.arange(65)))
        self.assertEqual(set(draw.tolist()), set(range(65)))


if __name__ == "__main__":
    unittest.main()
