#!/usr/bin/env python3
"""Synthetic tests for the bounded TeloHAEC extractor."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("telohaec_extract", HERE / "extract.py")
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_synthetic_h5ad(path: Path) -> np.ndarray:
    # Raw columns: A, duplicated DUP, and two panel genes P1/T1.
    rows = np.asarray(
        [
            [10, 2, 3, 99, 77],  # safe control c1
            [10, 2, 3, 88, 66],
            [20, 1, 4, 55, 44],  # non-targeting control c2
            [20, 1, 4, 33, 22],
            [40, 5, 5, 90, 90],  # T1 guide g1
            [40, 5, 5, 90, 90],
            [60, 7, 3, 80, 80],  # T1 guide g2
            [60, 7, 3, 80, 80],
            [900, 1, 1, 70, 70],  # T1 outlier guide g3
            [900, 1, 1, 70, 70],
            [30, 3, 2, 60, 60],  # compound target, excluded
            [30, 3, 2, 60, 60],
        ],
        dtype=np.float32,
    )
    matrix = sparse.csr_matrix(rows)
    barcodes = [
        "safe-targeting:CCCC:cell1-batch",
        "safe-targeting:CCCC:cell2-batch",
        "negative-control:DDDD:cell3-batch",
        "negative-control:DDDD:cell4-batch",
        "T1:AAAA:cell5-batch",
        "T1:AAAA:cell6-batch",
        "T1-TSS2:AAAT:cell7-batch",
        "T1-TSS2:AAAT:cell8-batch",
        "T1:AAAG:cell9-batch",
        "T1:AAAG:cell10-batch",
        "X-and-Y:GGGG:cell11-batch",
        "X-and-Y:GGGG:cell12-batch",
    ]
    features = ["A:ENSG1", "DUP:ENSG2", "DUP:ENSG3", "P1:ENSG4", "T1:ENSG5"]
    string_dtype = h5py.string_dtype("utf-8")
    with h5py.File(path, "w") as handle:
        x = handle.create_group("X")
        x.attrs["encoding-type"] = "csr_matrix"
        x.attrs["encoding-version"] = "0.1.0"
        x.attrs["shape"] = np.asarray(matrix.shape, dtype=np.int64)
        x.create_dataset("data", data=matrix.data, compression="gzip")
        x.create_dataset("indices", data=matrix.indices.astype(np.int32), compression="gzip")
        x.create_dataset("indptr", data=matrix.indptr.astype(np.int32), compression="gzip")
        obs = handle.create_group("obs")
        obs.create_dataset("barcodes", data=np.asarray(barcodes, dtype=object), dtype=string_dtype)
        var = handle.create_group("var")
        var.create_dataset("_index", data=np.asarray(features, dtype=object), dtype=string_dtype)
    return rows


def write_guides(path: Path) -> None:
    fields = ["guide_id", "spacer", "targeting", "type", "description"]
    rows = [
        ["T1__AAAA", "AAAA", "TRUE", "targeting", "T1"],
        ["T1_TSS2__AAAT", "AAAT", "TRUE", "targeting", "T1_TSS2"],
        ["T1__AAAG", "AAAG", "TRUE", "targeting", "T1"],
        ["safe_targeting__CCCC", "CCCC", "TRUE", "safe-targeting", "safe-targeting"],
        ["negative_control__DDDD", "DDDD", "FALSE", "non-targeting", "non-targeting"],
        ["X_and_Y__GGGG", "GGGG", "TRUE", "targeting", "X_and_Y"],
    ]
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(fields)
        writer.writerows(rows)


def write_panel(path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["target_gene"])
        writer.writerows([[f"P{i}"] for i in range(1, 300)])
        writer.writerow(["T1"])


class ExtractorTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory(dir=MODULE.DATA_ROOT, prefix=".synthetic-")
        self.root = Path(self.temp.name)
        self.h5ad = self.root / "synthetic.h5ad"
        self.guides = self.root / "guides.tsv.gz"
        self.panel = self.root / "panel.csv"
        self.rows = write_synthetic_h5ad(self.h5ad)
        write_guides(self.guides)
        write_panel(self.panel)

    def tearDown(self) -> None:
        self.temp.cleanup()

    def run_extract(self, name: str) -> Path:
        output = self.root / name
        MODULE.extract(
            self.h5ad,
            self.guides,
            self.panel,
            output,
            min_guide_cells=2,
            min_target_cells=4,
            row_block=3,
        )
        return output

    def test_target_mapping(self) -> None:
        self.assertEqual(MODULE.clean_target("ARHGEF15_TSS2"), ("ARHGEF15", ""))
        self.assertEqual(MODULE.clean_target("ARHGEF15-TSS2"), ("ARHGEF15", ""))
        self.assertEqual(MODULE.clean_target("Tm4sf1"), ("TM4SF1", ""))
        self.assertEqual(MODULE.clean_target("X_and_Y")[1], "compound_target")
        self.assertEqual(
            MODULE.clean_target("AL354709_1", "9p21_add_AL354709_1__AAAA")[1],
            "regulatory_element",
        )

    def test_output_scope_enforcement(self) -> None:
        with self.assertRaises(ValueError):
            MODULE.ensure_allowed_output(Path("/tmp/forbidden-telohaec-output"))

    def test_target_minimum_uses_only_eligible_guide_cells(self) -> None:
        guides = [
            MODULE.Guide("T1-a", "AAAA", "targeting", "T1", "T1", ""),
            MODULE.Guide("T1-b", "AAAT", "targeting", "T1", "T1", ""),
            MODULE.Guide("T1-c", "AAAG", "targeting", "T1", "T1", ""),
            MODULE.Guide("T2-a", "AATA", "targeting", "T2", "T2", ""),
            MODULE.Guide("T2-b", "AATT", "targeting", "T2", "T2", ""),
        ]
        counts = np.asarray([5, 4, 11, 5, 15], dtype=np.int64)
        mapping, eligible, total_cells, eligible_cells = MODULE._eligibility(
            guides,
            counts,
            min_guide_cells=5,
            min_target_cells=20,
        )

        self.assertEqual(total_cells, {"T1": 20, "T2": 20})
        self.assertEqual(eligible_cells, {"T1": 16, "T2": 20})
        self.assertEqual(eligible, [3, 4])
        np.testing.assert_array_equal(mapping, [-1, -1, -1, 0, 1])

    def test_extraction_math_and_invariants(self) -> None:
        output = self.run_extract("first")
        genes = (output / "genes.txt").read_text().splitlines()
        perts = (output / "perts.txt").read_text().splitlines()
        delta = np.load(output / "delta.npy")
        b_ctrl = np.load(output / "b_ctrl.npy")
        self.assertEqual(genes, ["A", "DUP", "P1", "T1"])
        self.assertEqual(perts, ["T1"])
        self.assertEqual(delta.shape, (1, 4))

        # Panel coordinates never enter normalization and remain exactly zero.
        self.assertTrue(np.array_equal(delta[:, [2, 3]], np.zeros((1, 2), dtype=np.float32)))
        self.assertTrue(np.array_equal(b_ctrl[[2, 3]], np.zeros(2, dtype=np.float32)))

        # Duplicate DUP columns are summed before per-guide normalization.
        c1 = np.asarray([20.0, 10.0])
        c2 = np.asarray([40.0, 10.0])
        expected_control = np.median(
            np.stack([MODULE.log_pseudobulk(c1), MODULE.log_pseudobulk(c2)]), axis=0
        )
        np.testing.assert_allclose(b_ctrl[:2], expected_control, rtol=1e-6)

        guide_profiles = np.stack(
            [
                MODULE.log_pseudobulk(np.asarray([80.0, 20.0])),
                MODULE.log_pseudobulk(np.asarray([120.0, 20.0])),
                MODULE.log_pseudobulk(np.asarray([1800.0, 4.0])),
            ]
        )
        expected_effect = np.median(guide_profiles, axis=0) - expected_control
        mean_effect = np.mean(guide_profiles, axis=0) - expected_control
        np.testing.assert_allclose(delta[0, :2], expected_effect, rtol=1e-6)
        self.assertFalse(np.allclose(delta[0, :2], mean_effect))

        report = json.loads((output / "extraction_report.json").read_text())
        self.assertEqual(report["schema"]["duplicate_feature_rows_collapsed"], 1)
        self.assertTrue(report["panel_coordinate_exclusion"]["delta_exact_zero_verified"])

    def test_deterministic_outputs(self) -> None:
        first = self.run_extract("deterministic-a")
        second = self.run_extract("deterministic-b")
        names = [
            "delta.npy",
            "b_ctrl.npy",
            "genes.txt",
            "perts.txt",
            "guide_qc.csv",
            "extraction_report.json",
        ]
        for name in names:
            self.assertEqual((first / name).read_bytes(), (second / name).read_bytes(), name)


if __name__ == "__main__":
    unittest.main(verbosity=2)
