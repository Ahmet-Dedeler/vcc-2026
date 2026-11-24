#!/usr/bin/env python3
"""Synthetic-only tests for Experiment 107 source-model input assembly."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
from scipy import sparse


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import model_input  # noqa: E402


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


class SyntheticModelInputFixture:
    """Small 16-channel fixture with the exact production artifact topology."""

    target_count = 3
    guide_count = 12
    response_count = 4
    exportable_count = 3
    vcc_mapped_response_count = 4
    rows_per_channel = 8

    def __init__(self, root: Path, *, forbidden_response_shift: int = 0) -> None:
        self.root = root
        self.work = root / "work"
        self.work.mkdir(parents=True)
        self.output = self.work / "source-model-input.npz"
        self.forbidden_response_shift = forbidden_response_shift
        self.source_verification = {
            "raw_archive_sha256": "a" * 64,
            "oracle_receipt_sha256": "b" * 64,
            "preflight_receipt_sha256": "c" * 64,
            "numpy_version": np.__version__,
            "scipy_version": model_input.source_artifacts.scipy.__version__,
        }
        self.guide_to_target = np.concatenate(
            (
                np.repeat(np.arange(self.target_count, dtype=np.int32), 3),
                np.full(3, -1, dtype=np.int32),
            )
        )
        self.crosswalks = self._crosswalks()
        self._write_channels()
        self._write_support_and_summary()

    def _crosswalks(self) -> SimpleNamespace:
        records = tuple(
            SimpleNamespace(
                symbol=f"R{index}",
                vcc_gene_index=index,
                is_panel_target=index == self.response_count - 1,
            )
            for index in range(self.response_count)
        )
        return SimpleNamespace(
            guide=SimpleNamespace(
                target_symbols=tuple(f"T{index}" for index in range(self.target_count)),
                guide_ids=tuple(f"G{index}" for index in range(self.guide_count)),
                guide_to_target=self.guide_to_target,
            ),
            response=SimpleNamespace(records=records),
        )

    def _target_design(self, channel: int) -> np.ndarray:
        rows = np.arange(self.rows_per_channel)
        return np.column_stack(
            (
                rows & 1,
                (rows >> 1) & 1,
                (rows >> 2) & 1,
            )
        ).astype(np.float64)

    def _write_channels(self) -> None:
        for channel in range(1, model_input.CHANNEL_COUNT + 1):
            directory = self.work / f"channel{channel:02d}"
            directory.mkdir()
            target_dense = self._target_design(channel)
            rows = np.arange(self.rows_per_channel)
            unmapped = np.column_stack(
                (
                    (rows + channel) % 3 == 0,
                    (rows * channel) % 5 < 2,
                    rows == channel % self.rows_per_channel,
                )
            ).astype(np.float64)
            guide_dense = np.column_stack(
                (np.repeat(target_dense, 3, axis=1), unmapped)
            )
            targets = sparse.csr_matrix(target_dense)
            guides = sparse.csr_matrix(guide_dense)
            library_size = np.full(
                self.rows_per_channel,
                1_000 + channel,
                dtype=np.int64,
            )
            response_counts = np.column_stack(
                (
                    5 + np.arange(self.rows_per_channel),
                    7 + target_dense[:, 0].astype(np.int64),
                    9 + target_dense[:, 1].astype(np.int64),
                    11
                    + target_dense[:, 2].astype(np.int64)
                    + self.forbidden_response_shift * (rows + 1),
                )
            )
            assignment_path = directory / "assignments.npz"
            np.savez_compressed(
                assignment_path,
                cell_indices=np.arange(self.rows_per_channel, dtype=np.int32),
                condition_codes=np.zeros(self.rows_per_channel, dtype=np.uint8),
                library_size=library_size,
                called_guide_burden=np.diff(guides.indptr).astype(np.uint16),
                non_targeting_guide_burden=unmapped[:, 0].astype(np.uint16),
                barcodes=np.asarray(
                    [
                        f"channel{channel:02d}-cell{row:02d}"
                        for row in range(self.rows_per_channel)
                    ],
                    dtype="U",
                ),
                target_indptr=targets.indptr.astype(np.int64),
                target_indices=targets.indices.astype(np.int32),
                called_guide_indptr=guides.indptr.astype(np.int64),
                called_guide_indices=guides.indices.astype(np.int32),
            )
            assignment_report_path = directory / "assignment-report.json"
            write_json(
                assignment_report_path,
                {
                    "status": "complete",
                    "output_sha256": model_input.sha256_file(assignment_path),
                    "assignment_stage_sha256": model_input.sha256_file(
                        model_input.HERE / "assignment_stage.py"
                    ),
                },
            )

            response_path = directory / "response-counts.npz"
            np.savez_compressed(
                response_path,
                cell_indices=np.arange(self.rows_per_channel, dtype=np.int32),
                response_counts=response_counts,
                library_size=library_size,
            )
            write_json(
                directory / "response-report.json",
                {
                    "status": "complete",
                    "assignment_report_sha256": model_input.sha256_file(
                        assignment_report_path
                    ),
                    "assignment_output_sha256": model_input.sha256_file(
                        assignment_path
                    ),
                    "output_sha256": model_input.sha256_file(response_path),
                    "assignment_stage_sha256": model_input.sha256_file(
                        model_input.HERE / "assignment_stage.py"
                    ),
                    "metadata_sha256": model_input.sha256_file(
                        model_input.HERE / "metadata.py"
                    ),
                    "source_artifacts_sha256": model_input.sha256_file(
                        model_input.HERE / "source_artifacts.py"
                    ),
                    "source_raw_archive_sha256": self.source_verification[
                        "raw_archive_sha256"
                    ],
                },
            )

    def _write_support_and_summary(self) -> None:
        target_counts = np.zeros((2, self.target_count), dtype=np.int64)
        guide_counts = np.zeros((2, self.guide_count), dtype=np.int64)
        target_channels = np.zeros(self.target_count, dtype=np.int16)
        for channel in range(1, model_input.CHANNEL_COUNT + 1):
            target_dense = self._target_design(channel)
            rows = np.arange(self.rows_per_channel)
            unmapped = np.column_stack(
                (
                    (rows + channel) % 3 == 0,
                    (rows * channel) % 5 < 2,
                    rows == channel % self.rows_per_channel,
                )
            ).astype(np.float64)
            guide_dense = np.column_stack(
                (np.repeat(target_dense, 3, axis=1), unmapped)
            )
            target_counts[0] += target_dense.sum(axis=0).astype(np.int64)
            guide_counts[0] += guide_dense.sum(axis=0).astype(np.int64)
            target_channels += target_dense.sum(axis=0) > 0
        distinct = np.full(self.target_count, 3, dtype=np.int32)
        guides_with_three = np.full(self.target_count, 3, dtype=np.int32)
        eligible = (
            (distinct >= model_input.modeling.MIN_DISTINCT_CALLED_GUIDES)
            & (
                guides_with_three
                >= model_input.modeling.MIN_GUIDES_WITH_THREE_UNTREATED_CELLS
            )
            & (target_counts[0] >= model_input.modeling.MIN_UNTREATED_CARRIER_CELLS)
            & (target_channels >= model_input.modeling.MIN_TRAINING_CHANNELS)
        )
        support_path = self.work / "target-support.npz"
        np.savez_compressed(
            support_path,
            guide_cells_before_qc=guide_counts,
            guide_cells_after_qc=guide_counts,
            target_cells_before_qc=target_counts,
            target_cells_after_qc=target_counts,
            distinct_called_guides=distinct,
            untreated_guides_with_at_least_three_cells=guides_with_three,
            untreated_support_channels=target_channels,
            eligible=eligible,
        )
        metadata_path = self.work / "metadata-crosswalk.npz"
        np.savez_compressed(metadata_path, sentinel=np.asarray([107]))
        metadata_report_path = self.work / "metadata-crosswalk-report.json"
        write_json(
            metadata_report_path,
            {"output_sha256": model_input.sha256_file(metadata_path)},
        )
        write_json(
            self.work / "assignment-summary.json",
            {
                "status": "pass",
                "gates": {
                    "retained_cells": True,
                    "untreated_fraction": True,
                    "median_called_guides": True,
                },
                "assignment_stage_sha256": model_input.sha256_file(
                    model_input.HERE / "assignment_stage.py"
                ),
                "source_artifacts_sha256": model_input.sha256_file(
                    model_input.HERE / "source_artifacts.py"
                ),
                "source_raw_archive_sha256": self.source_verification[
                    "raw_archive_sha256"
                ],
                "source_oracle_receipt_sha256": self.source_verification[
                    "oracle_receipt_sha256"
                ],
                "target_support_sha256": model_input.sha256_file(support_path),
                "metadata_crosswalk_report_sha256": model_input.sha256_file(
                    metadata_report_path
                ),
                "metadata_crosswalk_output_sha256": model_input.sha256_file(
                    metadata_path
                ),
                "channel_report_sha256": {
                    str(channel): model_input.sha256_file(
                        self.work / f"channel{channel:02d}" / "assignment-report.json"
                    )
                    for channel in range(1, model_input.CHANNEL_COUNT + 1)
                },
            },
        )

    def build(self) -> dict[str, object]:
        with (
            mock.patch.object(model_input, "WORK", self.work),
            mock.patch.object(
                model_input,
                "EXPECTED_EXPORTABLE_RESPONSE_COUNT",
                self.exportable_count,
            ),
            mock.patch.object(
                model_input,
                "EXPECTED_ORIGINAL_RESPONSE_COUNT",
                self.response_count,
            ),
            mock.patch.object(
                model_input,
                "EXPECTED_VCC_MAPPED_RESPONSE_COUNT",
                self.vcc_mapped_response_count,
            ),
            mock.patch.object(
                model_input.metadata,
                "load_pinned_crosswalks",
                return_value=self.crosswalks,
            ),
            mock.patch.object(
                model_input.source_artifacts,
                "verify_source_artifacts",
                return_value=self.source_verification,
            ),
        ):
            return model_input.build_model_input(self.output)


class ModelInputTests(unittest.TestCase):
    def test_csr_construction_preserves_empty_rows(self) -> None:
        matrix = model_input.csr_from_components(
            np.array([0, 2, 2, 3]),
            np.array([1, 4, 0]),
            rows=3,
            columns=5,
        )
        np.testing.assert_array_equal(
            matrix.toarray(),
            [
                [0, 1, 0, 0, 1],
                [0, 0, 0, 0, 0],
                [1, 0, 0, 0, 0],
            ],
        )

    def test_csr_construction_rejects_duplicates_and_unsorted_rows(self) -> None:
        for indices in (np.array([1, 1]), np.array([2, 1])):
            with self.assertRaisesRegex(model_input.ModelInputError, "canonical"):
                model_input.csr_from_components(
                    np.array([0, 2]),
                    indices,
                    rows=1,
                    columns=3,
                )

    def test_log_cpm_uses_full_library_size(self) -> None:
        counts = np.array([[10, 0], [5, 5]])
        library = np.array([100, 1_000])
        observed = model_input.log_cpm_response(counts, library)
        expected = np.log1p(
            np.array(
                [
                    [1_000.0, 0.0],
                    [50.0, 50.0],
                ]
            )
        )
        np.testing.assert_allclose(observed, expected, rtol=0, atol=0)

    def test_log_cpm_rejects_panel_counts_above_library(self) -> None:
        with self.assertRaisesRegex(model_input.ModelInputError, "exceed"):
            model_input.log_cpm_response(
                np.array([[6, 5]]),
                np.array([10]),
            )

    def test_assignment_summary_requires_exact_gates_and_current_reports(self) -> None:
        source_verification = {
            "raw_archive_sha256": "a" * 64,
            "oracle_receipt_sha256": "b" * 64,
        }
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            for channel in range(1, model_input.CHANNEL_COUNT + 1):
                channel_directory = work / f"channel{channel:02d}"
                channel_directory.mkdir()
                (channel_directory / "assignment-report.json").write_text(
                    f'{{"channel": {channel}}}\n',
                    encoding="utf-8",
                )
            (work / "target-support.npz").write_bytes(b"support")
            (work / "metadata-crosswalk.npz").write_bytes(b"metadata")
            metadata_output_sha256 = model_input.sha256_file(
                work / "metadata-crosswalk.npz"
            )
            (work / "metadata-crosswalk-report.json").write_text(
                ('{"output_sha256":"' + metadata_output_sha256 + '"}\n'),
                encoding="utf-8",
            )
            summary = {
                "status": "pass",
                "gates": {
                    "retained_cells": True,
                    "untreated_fraction": True,
                    "median_called_guides": True,
                },
                "assignment_stage_sha256": model_input.sha256_file(
                    model_input.HERE / "assignment_stage.py"
                ),
                "source_artifacts_sha256": model_input.sha256_file(
                    model_input.HERE / "source_artifacts.py"
                ),
                "source_raw_archive_sha256": "a" * 64,
                "source_oracle_receipt_sha256": "b" * 64,
                "target_support_sha256": model_input.sha256_file(
                    work / "target-support.npz"
                ),
                "metadata_crosswalk_report_sha256": model_input.sha256_file(
                    work / "metadata-crosswalk-report.json"
                ),
                "metadata_crosswalk_output_sha256": metadata_output_sha256,
                "channel_report_sha256": {
                    str(channel): model_input.sha256_file(
                        work / f"channel{channel:02d}" / "assignment-report.json"
                    )
                    for channel in range(1, model_input.CHANNEL_COUNT + 1)
                },
            }
            with mock.patch.object(model_input, "WORK", work):
                observed = model_input._validate_assignment_summary(
                    summary,
                    source_verification=source_verification,
                )
                self.assertEqual(
                    set(observed),
                    {
                        str(channel)
                        for channel in range(1, model_input.CHANNEL_COUNT + 1)
                    },
                )
                (work / "channel03" / "assignment-report.json").write_text(
                    '{"channel": 3, "changed": true}\n',
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(
                    model_input.ModelInputError,
                    "changed after assay summarization",
                ):
                    model_input._validate_assignment_summary(
                        summary,
                        source_verification=source_verification,
                    )

    def test_assignment_summary_rejects_missing_or_extra_gates(self) -> None:
        source_verification = {
            "raw_archive_sha256": "a" * 64,
            "oracle_receipt_sha256": "b" * 64,
        }
        for gates in ({}, {"retained_cells": True}, {"extra": True}):
            with self.subTest(gates=gates):
                with self.assertRaisesRegex(
                    model_input.ModelInputError,
                    "assay-compatibility",
                ):
                    model_input._validate_assignment_summary(
                        {"status": "pass", "gates": gates},
                        source_verification=source_verification,
                    )

    def test_build_model_input_production_path_carries_exact_eligibility(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as directory:
            fixture = SyntheticModelInputFixture(Path(directory))
            report = fixture.build()
            self.assertEqual(report["status"], "complete")
            self.assertEqual(
                report["source_preflight_receipt_sha256"],
                fixture.source_verification["preflight_receipt_sha256"],
            )
            with np.load(fixture.output, allow_pickle=False) as values:
                self.assertIn("precomputed_eligible", values.files)
                self.assertTrue(np.all(values["precomputed_eligible"]))
                self.assertEqual(
                    values["response"].shape,
                    (
                        model_input.CHANNEL_COUNT * fixture.rows_per_channel,
                        fixture.exportable_count,
                    ),
                )
                np.testing.assert_array_equal(
                    values["response_symbols"],
                    ["R0", "R1", "R2"],
                )
                np.testing.assert_array_equal(
                    values["response_exportable"],
                    np.ones(fixture.exportable_count, dtype=bool),
                )
                np.testing.assert_array_equal(
                    values["original_response_symbols"],
                    ["R0", "R1", "R2", "R3"],
                )
                np.testing.assert_array_equal(
                    values["original_response_exportable"],
                    [True, True, True, False],
                )
                np.testing.assert_array_equal(
                    values["response_original_indices"],
                    [0, 1, 2],
                )
                self.assertEqual(
                    values["precomputed_distinct_called_guides"].shape,
                    (fixture.target_count,),
                )
            audit_path, audit_report_path = model_input.design_audit_paths(
                fixture.output
            )
            self.assertTrue(audit_path.is_file())
            self.assertTrue(audit_report_path.is_file())
            audit_report = json.loads(audit_report_path.read_text(encoding="utf-8"))
            self.assertEqual(audit_report["status"], "complete_before_response_fit")
            self.assertIs(audit_report["response_values_used"], False)
            self.assertEqual(
                report["design_audit_output_sha256"],
                model_input.sha256_file(audit_path),
            )

    def test_design_audit_rejects_tampered_nuisance_binding(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as directory:
            fixture = SyntheticModelInputFixture(Path(directory))
            fixture.build()
            with np.load(fixture.output, allow_pickle=False) as values:
                arrays = {name: values[name] for name in values.files}
            target_design = model_input.csr_from_components(
                arrays["target_indptr"],
                arrays["target_indices"],
                rows=arrays["response"].shape[0],
                columns=fixture.target_count,
            )
            guide_presence = model_input.csr_from_components(
                arrays["guide_indptr"],
                arrays["guide_indices"],
                rows=arrays["response"].shape[0],
                columns=fixture.guide_count,
            )
            audit_path, _ = model_input.design_audit_paths(fixture.output)
            with self.assertRaisesRegex(
                model_input.ModelInputError,
                "design_binding_sha256",
            ):
                with mock.patch.object(
                    model_input,
                    "EXPECTED_EXPORTABLE_RESPONSE_COUNT",
                    fixture.exportable_count,
                ):
                    model_input.load_production_design_audit(
                        audit_path,
                        target_design=target_design,
                        nuisance=np.asarray(arrays["nuisance"], dtype=np.float64) + 1.0,
                        channels=arrays["channels"],
                        guide_evidence=model_input.modeling.GuideEvidence(
                            guide_presence=guide_presence,
                            guide_target_indices=arrays["guide_to_target"],
                        ),
                    )

    def test_build_model_input_rejects_precomputed_eligibility_drift(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as directory:
            fixture = SyntheticModelInputFixture(Path(directory))
            support_path = fixture.work / "target-support.npz"
            with np.load(support_path, allow_pickle=False) as values:
                arrays = {name: values[name] for name in values.files}
            arrays["eligible"] = arrays["eligible"].copy()
            arrays["eligible"][0] = False
            np.savez_compressed(support_path, **arrays)
            summary_path = fixture.work / "assignment-summary.json"
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            summary["target_support_sha256"] = model_input.sha256_file(support_path)
            write_json(summary_path, summary)
            with self.assertRaisesRegex(
                model_input.ModelInputError,
                "eligibility is inconsistent",
            ):
                fixture.build()


if __name__ == "__main__":
    unittest.main()
