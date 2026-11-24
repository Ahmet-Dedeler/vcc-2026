"""Small tests for experiment 049; no real-data load and no cell-eval2."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from run import (  # noqa: E402
    PRIMARY_METRICS,
    arm_summary,
    emit_template_fixed_margins,
    macro_summary,
    scorer_pseudobulk,
    target_seed,
)


class MultiTargetSummaryTest(unittest.TestCase):
    def test_arm_summary_extracts_requested_metrics(self) -> None:
        got = arm_summary(
            library_ks=0.2,
            detection={"gene_detection_rate_mae_vs_real_target": 0.1},
            cpm={"lfc_nmae_proxy_vs_real": 1.3},
            de={
                "bh05_jaccard_vs_real": 0.04,
                "bh05_direction_accuracy_in_intersection": 0.6,
                "bh05_n_sig": 100,
                "real_bh05_n_sig": 50,
                "top_k": {"500": {"set_overlap_fraction": 0.12}},
            },
        )
        self.assertEqual(set(PRIMARY_METRICS), set(got) - {
            "predicted_bh05_n_sig",
            "real_bh05_n_sig",
        })
        self.assertEqual(got["top500_set_overlap_fraction"], 0.12)

    def test_macro_summary_uses_equal_target_weights_and_directions(self) -> None:
        per_target = {}
        for i, target in enumerate(("A", "B", "C"), start=1):
            quasi = {
                metric: float(i)
                for metric in PRIMARY_METRICS
            }
            template = dict(quasi)
            for metric, direction in PRIMARY_METRICS.items():
                template[metric] += 1.0 if direction == "higher" else -1.0
            per_target[target] = {
                "arms": {"quasi": quasi, "template": template}
            }

        got = macro_summary(per_target)
        for metric in PRIMARY_METRICS:
            self.assertEqual(got["arms"]["quasi"][metric], 2.0)
            self.assertEqual(
                got["comparison"][metric]["favorable_delta_for_template"], 1.0
            )
            self.assertEqual(
                got["comparison"][metric]["template_target_wins"], 3
            )

    def test_template_integration_preserves_scorer_pseudobulk(self) -> None:
        template = sp.csr_matrix(
            np.array(
                [
                    [4, 0, 2, 0],
                    [0, 5, 0, 1],
                    [2, 0, 3, 4],
                    [0, 3, 1, 2],
                ],
                dtype=np.int64,
            )
        )
        columns = np.array([8, 10, 7, 6], dtype=np.int64)
        rows = np.array([7, 7, 9, 8], dtype=np.int64)
        emitted, _ = emit_template_fixed_margins(
            template,
            column_totals=columns,
            row_margins=rows,
            seed=49,
        )
        np.testing.assert_array_equal(
            np.asarray(emitted.sum(axis=0)).ravel(), columns
        )
        np.testing.assert_array_equal(
            scorer_pseudobulk(emitted),
            np.log1p(50_000.0 * columns / columns.sum()),
        )

    def test_med25_reuses_experiment_043_seed(self) -> None:
        self.assertEqual(target_seed(4300, "MED25", 1), 4300)
        self.assertEqual(target_seed(4300, "MED15", 0), 4400)


if __name__ == "__main__":
    unittest.main()
