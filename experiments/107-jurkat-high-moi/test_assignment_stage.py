#!/usr/bin/env python3
"""Synthetic-only tests for Experiment 107's streaming assignment stage."""

from __future__ import annotations

import io
import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import assignment_stage  # noqa: E402


def matrix_stream(
    rows: int,
    columns: int,
    entries: list[tuple[int, int, int]],
) -> io.BytesIO:
    lines = [
        "%%MatrixMarket matrix coordinate integer general",
        "% synthetic",
        f"{rows} {columns} {len(entries)}",
        *(f"{row + 1} {column + 1} {value}" for row, column, value in entries),
    ]
    return io.BytesIO(("\n".join(lines) + "\n").encode())


class AssignmentStageTests(unittest.TestCase):
    def test_dense_extraction_keeps_selected_column_order(self) -> None:
        entries = [
            (0, 0, 2),
            (2, 0, 1),
            (1, 2, 5),
            (0, 4, 3),
            (2, 4, 7),
        ]
        observed = assignment_stage.extract_selected_dense(
            matrix_stream(3, 5, entries),
            np.array([0, 3, 4]),
            expected_rows=3,
            expected_columns=5,
        )
        np.testing.assert_array_equal(
            observed,
            np.array(
                [
                    [2, 0, 1],
                    [0, 0, 0],
                    [3, 0, 7],
                ]
            ),
        )

    def test_sparse_batches_include_empty_selected_cells(self) -> None:
        entries = [
            (0, 0, 2),
            (2, 0, 1),
            (1, 2, 5),
            (0, 4, 3),
            (2, 4, 7),
        ]
        batches = list(
            assignment_stage.iter_selected_sparse_batches(
                matrix_stream(3, 5, entries),
                np.array([0, 1, 2, 4]),
                expected_rows=3,
                expected_columns=5,
                batch_cells=2,
            )
        )
        self.assertEqual(len(batches), 2)
        np.testing.assert_array_equal(batches[0].cell_indices, [0, 1])
        np.testing.assert_array_equal(batches[0].indptr, [0, 2, 2])
        np.testing.assert_array_equal(batches[0].guide_indices, [0, 2])
        np.testing.assert_array_equal(batches[0].observed_umi, [2, 1])
        np.testing.assert_array_equal(batches[1].cell_indices, [2, 4])
        np.testing.assert_array_equal(batches[1].indptr, [0, 1, 3])
        np.testing.assert_array_equal(batches[1].guide_indices, [1, 0, 2])
        np.testing.assert_array_equal(batches[1].observed_umi, [5, 3, 7])

    def test_streams_reject_non_column_major_coordinates(self) -> None:
        entries = [(0, 1, 1), (1, 0, 1)]
        with self.assertRaisesRegex(
            assignment_stage.AssignmentStageError,
            "column-major",
        ):
            assignment_stage.extract_selected_dense(
                matrix_stream(2, 2, entries),
                np.array([0]),
                expected_rows=2,
                expected_columns=2,
            )

    def test_expression_qc_uses_frozen_inclusive_thresholds(self) -> None:
        totals = np.array([1_000, 1_100, 1_200, 1, 1_300])
        detected = np.array([10, 9, 10, 10, 10])
        mitochondrial = np.array([300, 0, 361, 0, 0])
        keep, report = assignment_stage.expression_qc_mask(
            totals,
            detected,
            mitochondrial,
        )
        np.testing.assert_array_equal(keep, [True, False, False, False, True])
        self.assertEqual(report["failed_detected_genes"], 1)
        self.assertEqual(report["failed_mitochondrial_fraction"], 1)
        self.assertEqual(report["failed_log_library_mad"], 1)

    def test_barcode_extraction_validates_full_axis(self) -> None:
        stream = io.BytesIO(b"A\nB\nC\nD\n")
        observed = assignment_stage.extract_selected_barcodes(
            stream,
            np.array([0, 2]),
            expected_columns=4,
        )
        np.testing.assert_array_equal(observed, ["A", "C"])
        with self.assertRaisesRegex(
            assignment_stage.AssignmentStageError,
            "row count",
        ):
            assignment_stage.extract_selected_barcodes(
                io.BytesIO(b"A\nB\n"),
                np.array([0]),
                expected_columns=4,
            )

    def test_response_extraction_sums_duplicate_symbol_rows(self) -> None:
        records = tuple(
            assignment_stage.metadata.ResponseGeneCrosswalk(
                response_index=index,
                symbol=f"G{index}",
                entrez_id=index + 1,
                group="synthetic",
                source_feature_indices=(0, 2) if index == 0 else (),
                vcc_gene_indices=(),
                is_panel_target=False,
            )
            for index in range(assignment_stage.metadata.EXPECTED_RESPONSE_GENES)
        )
        counts = assignment_stage.extract_response_counts(
            matrix_stream(
                3,
                2,
                [(0, 0, 2), (2, 0, 3), (1, 1, 11)],
            ),
            np.array([0]),
            records,
            expected_rows=3,
            expected_columns=2,
        )
        self.assertEqual(counts.shape, (1, 374))
        self.assertEqual(counts[0, 0], 5)
        self.assertEqual(int(counts.sum()), 5)

    def test_selected_columns_must_be_strict_and_in_bounds(self) -> None:
        for selected in (
            np.array([1, 1]),
            np.array([2, 1]),
            np.array([-1]),
            np.array([5]),
        ):
            with self.assertRaises(assignment_stage.AssignmentStageError):
                assignment_stage.extract_selected_dense(
                    matrix_stream(1, 5, []),
                    selected,
                    expected_rows=1,
                    expected_columns=5,
                )

    def test_target_support_thresholds_are_exact(self) -> None:
        guide_counts = np.zeros((2, 8), dtype=np.int64)
        mapping = np.array([0, 0, 0, 1, 1, 1, 2, -1])
        guide_counts[0, :3] = [3, 3, 1]
        guide_counts[1, 2] = 1
        guide_counts[0, 3:6] = [3, 3, 1]
        guide_counts[0, 6] = 100
        target_counts = np.array(
            [
                [20, 19, 100],
                [5, 5, 0],
            ]
        )
        channels = np.array([4, 4, 16])
        observed = assignment_stage.derive_target_support(
            guide_counts,
            target_counts,
            channels,
            mapping,
        )
        np.testing.assert_array_equal(
            observed["distinct_called_guides"],
            [3, 3, 1],
        )
        np.testing.assert_array_equal(
            observed["untreated_guides_with_at_least_three_cells"],
            [2, 2, 1],
        )
        np.testing.assert_array_equal(observed["eligible"], [True, False, False])

    def test_activated_only_guides_do_not_satisfy_untreated_support(self) -> None:
        guide_counts = np.zeros((2, 4), dtype=np.int64)
        guide_counts[0, :2] = [4, 4]
        guide_counts[1, 2] = 8
        observed = assignment_stage.derive_target_support(
            guide_counts,
            np.array([[24], [8]], dtype=np.int64),
            np.array([4], dtype=np.int16),
            np.array([0, 0, 0, -1], dtype=np.int64),
        )
        np.testing.assert_array_equal(observed["distinct_called_guides"], [2])
        np.testing.assert_array_equal(observed["eligible"], [False])

    def test_assignment_validation_recomputes_burdens_and_targets(self) -> None:
        arrays = {
            "cell_indices": np.array([2, 8], dtype=np.int32),
            "barcodes": np.array(["A", "B"]),
            "hash_label_indices": np.array([0, 1], dtype=np.int8),
            "condition_codes": np.array([0, 1], dtype=np.uint8),
            "library_size": np.array([100, 200]),
            "detected_genes": np.array([20, 30]),
            "mitochondrial_counts": np.array([0, 10]),
            "total_guide_umi": np.array([4, 5]),
            "called_guide_burden": np.array([2, 1]),
            "non_targeting_guide_burden": np.array([1, 0]),
            "called_guide_indptr": np.array([0, 2, 3]),
            "called_guide_indices": np.array([0, 2, 1]),
            "target_indptr": np.array([0, 1, 2]),
            "target_indices": np.array([1, 0]),
        }
        mapping = np.full(assignment_stage.EXPECTED_GUIDE_ROWS, -1)
        mapping[1] = 0
        mapping[2] = 1
        ntc = np.zeros(assignment_stage.EXPECTED_GUIDE_ROWS, dtype=bool)
        ntc[0] = True
        assignment_stage.validate_assignment_arrays(
            arrays,
            guide_to_target=mapping,
            non_targeting_guides=ntc,
            hash_condition_codes=np.array([0, 1]),
        )
        arrays["target_indices"] = np.array([0, 0])
        with self.assertRaisesRegex(
            assignment_stage.AssignmentStageError,
            "target CSR",
        ):
            assignment_stage.validate_assignment_arrays(
                arrays,
                guide_to_target=mapping,
                non_targeting_guides=ntc,
                hash_condition_codes=np.array([0, 1]),
            )

    def test_metadata_crosswalk_report_is_exact_and_reloadable(self) -> None:
        crosswalks = assignment_stage.metadata.load_pinned_crosswalks()
        old_work = assignment_stage.WORK
        try:
            with self.subTest("write"):
                import tempfile

                with tempfile.TemporaryDirectory() as directory:
                    assignment_stage.WORK = Path(directory)
                    report = assignment_stage.write_metadata_crosswalk_report(
                        crosswalks
                    )
                    self.assertEqual(report["guides"], 83_401)
                    self.assertEqual(report["source_targets"], 18_465)
                    self.assertEqual(report["response_genes"], 374)
                    self.assertEqual(report["pds_visible_response_genes"], 349)
                    path = Path(directory) / "metadata-crosswalk.npz"
                    self.assertEqual(
                        report["output_sha256"],
                        assignment_stage.sha256_file(path),
                    )
                    with np.load(path, allow_pickle=False) as values:
                        self.assertEqual(values["guide_to_target"].shape, (83_401,))
                        self.assertEqual(
                            values["response_source_indptr"].shape,
                            (375,),
                        )
                        self.assertEqual(
                            np.count_nonzero(
                                (values["response_vcc_indices"] >= 0)
                                & ~values["response_is_panel_target"]
                            ),
                            349,
                        )
        finally:
            assignment_stage.WORK = old_work


if __name__ == "__main__":
    unittest.main()
