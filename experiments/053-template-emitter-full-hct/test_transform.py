"""Synthetic integration tests for the full-HCT streaming transformer."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from compare_scores import collect
from transform import read_rows, transform_file

GENES = ["G1", "G2", "G3", "G4", "G5", "G6"]
PERTURBATIONS = ("P1", "P2")
CONTROL = "non-targeting"


def fixture_counts() -> tuple[sp.csr_matrix, sp.csr_matrix, list[str]]:
    """Return small prediction/real matrices with a valid contiguous layout."""
    prediction = np.array(
        [
            [8, 0, 1, 0, 2, 0],
            [7, 1, 0, 0, 3, 0],
            [9, 0, 2, 1, 0, 0],
            [6, 2, 0, 0, 1, 1],
            [0, 7, 2, 0, 0, 1],
            [1, 8, 0, 1, 0, 0],
            [0, 6, 3, 0, 1, 0],
            [2, 9, 0, 0, 0, 1],
            [4, 1, 0, 2, 0, 1],
            [2, 3, 1, 0, 1, 1],
            [6, 0, 2, 0, 1, 0],
            [0, 5, 0, 2, 0, 2],
            [3, 1, 1, 1, 1, 1],
            [1, 4, 0, 1, 2, 0],
            [5, 0, 1, 2, 0, 1],
            [0, 3, 2, 1, 1, 2],
        ],
        dtype=np.int32,
    )
    real = prediction.copy()
    real[:4, [0, 3]] += np.array([[-2, 2]])
    real[4:8, [1, 4]] += np.array([[-2, 2]])
    labels = ["P1"] * 4 + ["P2"] * 4 + [CONTROL] * 8
    return sp.csr_matrix(prediction), sp.csr_matrix(real), labels


def write_fixture(
    prediction_path: Path,
    real_path: Path,
) -> None:
    prediction, real, labels = fixture_counts()
    obs = pd.DataFrame(
        {
            "target_gene": pd.Categorical(labels),
            "context": pd.Categorical(["HCT116"] * len(labels)),
        },
        index=[f"cell-{i}" for i in range(len(labels))],
    )
    var = pd.DataFrame(
        {"feature_type": ["Gene Expression"] * len(GENES)},
        index=GENES,
    )
    pred = ad.AnnData(X=prediction, obs=obs.copy(), var=var.copy())
    pred.uns["fixture"] = "prediction"
    truth = ad.AnnData(X=real, obs=obs.copy(), var=var.copy())
    truth.uns["fixture"] = "real"
    pred.write_h5ad(prediction_path, compression="lzf")
    truth.write_h5ad(real_path, compression="lzf")


class StreamingTransformerTest(unittest.TestCase):
    def test_full_synthetic_transform_preserves_contract_and_pseudobulk(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prediction_path = root / "prediction.h5ad"
            real_path = root / "real.h5ad"
            output_path = root / "template.h5ad"
            report_path = root / "report.json"
            write_fixture(prediction_path, real_path)

            report = transform_file(
                prediction_path,
                real_path,
                output_path,
                report_path=report_path,
                cells_per_perturbation=4,
                template_cells=4,
                chunk_rows=2,
                chunk_nnz=16,
                max_rss_gib=3.5,
            )

            self.assertTrue(output_path.is_file())
            self.assertTrue(report_path.is_file())
            self.assertEqual(report["perturbations_transformed"], 2)
            self.assertTrue(report["all_in_memory_pseudobulk_exact"])
            self.assertTrue(report["post_write_verification"]["all_pseudobulk_exact"])
            self.assertTrue(report["post_write_verification"]["controls_exact"])

            source = ad.read_h5ad(prediction_path, backed="r")
            output = ad.read_h5ad(output_path, backed="r")
            try:
                self.assertEqual(source.shape, output.shape)
                self.assertTrue(source.obs.equals(output.obs))
                self.assertTrue(source.var.equals(output.var))
                self.assertEqual(output.uns["fixture"], "prediction")
                labels = source.obs["target_gene"].astype(str).to_numpy()
                for target in PERTURBATIONS:
                    rows = np.flatnonzero(labels == target)
                    before = read_rows(source, rows)
                    after = read_rows(output, rows)
                    np.testing.assert_array_equal(
                        np.asarray(before.sum(axis=0)).ravel(),
                        np.asarray(after.sum(axis=0)).ravel(),
                    )
                control_rows = np.flatnonzero(labels == CONTROL)
                before = read_rows(source, control_rows)
                after = read_rows(output, control_rows)
                self.assertEqual((before != after).nnz, 0)
            finally:
                source.file.close()
                output.file.close()

    def test_dry_run_exercises_emitter_without_writing_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prediction_path = root / "prediction.h5ad"
            real_path = root / "real.h5ad"
            output_path = root / "must-not-exist.h5ad"
            write_fixture(prediction_path, real_path)

            report = transform_file(
                prediction_path,
                real_path,
                output_path,
                cells_per_perturbation=4,
                template_cells=4,
                chunk_rows=2,
                max_rss_gib=3.5,
                dry_run=True,
            )

            self.assertEqual(report["mode"], "dry-run")
            self.assertTrue(report["exact_pseudobulk"])
            self.assertFalse(output_path.exists())

    def test_score_collector_compares_paired_mean_rows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tag = root / "g0.65"
            template_run = tag / "template-run"
            template_run.mkdir(parents=True)
            fields = [
                "statistic",
                "pds_cosine",
                "de_wilcoxon_lfc_nmae",
            ]
            for path, values in (
                (
                    tag / "quasi-agg_results.csv",
                    ["mean", "0.70", "1.10"],
                ),
                (
                    template_run / "agg_results.csv",
                    ["mean", "0.71", "1.05"],
                ),
            ):
                with path.open("w", newline="") as handle:
                    writer = csv.writer(handle)
                    writer.writerow(fields)
                    writer.writerow(values)

            result = collect(root)
            comparison = result["runs"]["g0.65"]
            self.assertAlmostEqual(
                comparison["pds_cosine"]["favorable_delta_for_template"],
                0.01,
            )
            self.assertAlmostEqual(
                comparison["de_wilcoxon_lfc_nmae"]["favorable_delta_for_template"],
                0.05,
            )


if __name__ == "__main__":
    unittest.main()
