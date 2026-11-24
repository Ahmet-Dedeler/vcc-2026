from __future__ import annotations

import hashlib
import gzip
import json
import sys
import tempfile
import unittest
from array import array
from dataclasses import replace
from pathlib import Path
from typing import Any, Callable

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import source_artifacts  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def bh(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values, kind="stable")
    ranked = p_values[order] * p_values.size / np.arange(1, p_values.size + 1)
    adjusted_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(p_values)
    adjusted[order] = np.minimum(adjusted_sorted, 1.0)
    return adjusted


class SyntheticFixture:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.paths = source_artifacts.SourceArtifactPaths.from_root(root)
        self.channels = (1, 2)
        self.accessions = ("GSM0000001", "GSM0000002")
        self.columns = 6
        self.rows = {
            "transcriptome": 4,
            "guides": 3,
            "labels": 2,
        }
        self.paths.archive.parent.mkdir(parents=True)
        self.paths.work.mkdir(parents=True)
        (root / "audit" / "support").mkdir(parents=True)

        self.paths.archive.write_bytes(b"synthetic raw archive bytes")
        self.paths.file_inventory.write_text(
            "#Archive/File\tName\tTime\tSize\tType\n",
            encoding="utf-8",
        )
        self.paths.preprocess_code.write_text(
            "def preprocess():\n    return 107\n",
            encoding="utf-8",
        )
        self.paths.r_oracle_code.write_text(
            "stopifnot(TRUE)\n",
            encoding="utf-8",
        )
        self.paths.source_stage_code.write_text(
            "def source_stage():\n    return 107\n",
            encoding="utf-8",
        )
        self.paths.preflight_code.write_text(
            "def source_preflight():\n    return 107\n",
            encoding="utf-8",
        )
        for filename in (
            "guide_calling.py",
            "demultiplex.py",
            "test_preprocess.py",
            "test_source_stage.py",
            "test_guide_calling.py",
            "test_demultiplex.py",
        ):
            (root / filename).write_text(
                f"# synthetic {filename}\n",
                encoding="utf-8",
            )

        self.manifest = self._manifest()
        write_json(self.paths.manifest, self.manifest)
        self.provenance = {
            "status": "preregistered_before_acquisition",
            "source": {
                "raw_archive_size_bytes": self.paths.archive.stat().st_size,
            },
        }
        write_json(self.paths.provenance, self.provenance)
        self.protocol = self._protocol()
        write_json(self.paths.protocol, self.protocol)

        self.contract = replace(
            source_artifacts.PRODUCTION_CONTRACT,
            protocol_sha256=sha256(self.paths.protocol),
            manifest_sha256=sha256(self.paths.manifest),
            provenance_sha256=sha256(self.paths.provenance),
            raw_archive_sha256=sha256(self.paths.archive),
            raw_archive_size_bytes=self.paths.archive.stat().st_size,
            file_inventory_sha256=sha256(self.paths.file_inventory),
            channels=self.channels,
            accessions=self.accessions,
            barcode_columns=self.columns,
            feature_rows=tuple(self.rows.items()),
        )
        self._write_oracle()
        self._write_preflight()
        self._write_audits()
        for channel, accession in zip(
            self.channels,
            self.accessions,
            strict=True,
        ):
            self._write_channel(channel, accession)

    def _manifest(self) -> dict[str, Any]:
        return {
            "schema_version": 1,
            "dataset": "GSE249595",
            "assets": {
                "raw_archive": {
                    "filename": self.paths.archive.name,
                    "size_bytes": self.paths.archive.stat().st_size,
                    "sha256": None,
                },
                "file_inventory": {
                    "filename": self.paths.file_inventory.name,
                    "size_bytes": self.paths.file_inventory.stat().st_size,
                    "sha256": sha256(self.paths.file_inventory),
                },
            },
            "archive_contract": {
                "expected_member_count": len(self.channels) * 9,
                "channels": [
                    {"channel": channel, "accession": accession}
                    for channel, accession in zip(
                        self.channels,
                        self.accessions,
                        strict=True,
                    )
                ],
                "modalities": {
                    modality: {
                        "feature_rows": rows,
                        "barcode_columns": self.columns,
                    }
                    for modality, rows in self.rows.items()
                },
                "roles": {
                    "barcode": ".tsv.gz",
                    "features": ".tsv.gz",
                    "matrix": ".mtx.gz",
                },
            },
            "safety": {
                "max_assignment_nnz_total": 1_000,
                "assignment_coordinate_order": "row_or_column_major_strict",
                "never_scan_transcriptome_coordinates": True,
            },
        }

    def _protocol(self) -> dict[str, Any]:
        semantic = {
            "repository": "https://example.invalid/DropletUtils",
            "commit": "d" * 40,
            "package_version": "1.31.2",
            "empty_drops_source_sha256": "1" * 64,
            "barcode_ranks_source_sha256": "2" * 64,
            "ambient_profile_source_sha256": "3" * 64,
            "monte_carlo_source_sha256": "4" * 64,
            "good_turing_reference": {
                "repository": "https://example.invalid/edgeR",
                "commit": "e" * 40,
                "package_version": "4.99.0",
                "r_source_sha256": "5" * 64,
                "c_source_sha256": "6" * 64,
            },
        }
        return {
            "schema_version": 1,
            "experiment": "107-jurkat-high-moi",
            "status": (
                "acquisition_complete_and_protocol_frozen_before_"
                "transcriptome_coordinate_values_are_parsed"
            ),
            "cell_and_assignment_freeze": {
                "transcriptome_cell_call": {
                    "rng_contract": {
                        "numpy_version": np.__version__,
                    },
                },
                "guide_assignment": {
                    "scipy_version": source_artifacts.scipy.__version__,
                },
            },
            "immutable_inputs": {
                "acquisition_manifest": {
                    "path": "ACQUISITION_MANIFEST.json",
                    "sha256": sha256(self.paths.manifest),
                },
                "provenance": {
                    "path": "PROVENANCE.json",
                    "sha256": sha256(self.paths.provenance),
                },
                "emptydrops_semantic_reference": semantic,
            },
        }

    def _write_oracle(self) -> None:
        semantic = self.protocol["immutable_inputs"]["emptydrops_semantic_reference"]
        good_turing = semantic["good_turing_reference"]
        write_json(
            self.paths.oracle_receipt,
            {
                "status": "pass",
                "checks": list(self.contract.required_oracle_checks),
                "DropletUtils_version": semantic["package_version"],
                "DropletUtils_commit": semantic["commit"],
                "DropletUtils_empty_drops_sha256": semantic[
                    "empty_drops_source_sha256"
                ],
                "DropletUtils_barcode_ranks_sha256": semantic[
                    "barcode_ranks_source_sha256"
                ],
                "DropletUtils_ambient_profile_sha256": semantic[
                    "ambient_profile_source_sha256"
                ],
                "edgeR_version": good_turing["package_version"],
                "edgeR_commit": good_turing["commit"],
                "edgeR_good_turing_r_sha256": good_turing["r_source_sha256"],
                "edgeR_good_turing_c_sha256": good_turing["c_source_sha256"],
                "python_port_sha256": sha256(self.paths.preprocess_code),
                "r_oracle_sha256": sha256(self.paths.r_oracle_code),
                "synthetic_shape": [4, 8],
                "synthetic_nnz": 12,
            },
        )

    def _write_preflight(self) -> None:
        write_json(
            self.paths.preflight_receipt,
            {
                "schema_version": 1,
                "status": "pass",
                "numpy_version": np.__version__,
                "scipy_version": source_artifacts.scipy.__version__,
                "python_version": sys.version.split()[0],
                "protocol_sha256": sha256(self.paths.protocol),
                "preflight_implementation_sha256": sha256(
                    self.paths.preflight_code
                ),
                "implementation_sha256": {
                    filename: sha256(self.root / filename)
                    for filename in (
                        "source_stage.py",
                        "preprocess.py",
                        "guide_calling.py",
                        "demultiplex.py",
                    )
                },
                "test_module_sha256": {
                    filename: sha256(self.root / filename)
                    for filename in (
                        "test_preprocess.py",
                        "test_source_stage.py",
                        "test_guide_calling.py",
                        "test_demultiplex.py",
                    )
                },
                "tests_run": 39,
                "failures": 0,
                "errors": 0,
                "skipped": 0,
                "transcript_sha256": "7" * 64,
                "held_outcomes_opened": False,
            },
        )

    def _binding(self) -> dict[str, Any]:
        return {
            "manifest_sha256": sha256(self.paths.manifest),
            "provenance_sha256": sha256(self.paths.provenance),
            "raw_archive_sha256": sha256(self.paths.archive),
            "file_inventory_sha256": sha256(self.paths.file_inventory),
            "raw_hash_recomputed_for_this_run": True,
        }

    def _inventory(self) -> dict[str, Any]:
        return {
            "status": "pass",
            "archive": {
                "filename": self.paths.archive.name,
                "size_bytes": self.paths.archive.stat().st_size,
            },
            "file_inventory_sha256": sha256(self.paths.file_inventory),
            "member_count": len(self.channels) * 9,
            "member_total_compressed_bytes": 123,
            "channels": [
                {
                    "channel": channel,
                    "accession": accession,
                    "member_count": 9,
                }
                for channel, accession in zip(
                    self.channels,
                    self.accessions,
                    strict=True,
                )
            ],
        }

    def _write_audits(self) -> None:
        feature_hashes = {
            modality: hashlib.sha256(modality.encode()).hexdigest()
            for modality in self.rows
        }
        barcode_hash = hashlib.sha256(b"barcodes").hexdigest()
        support_records: list[dict[str, Any]] = []
        schema_channels: list[dict[str, Any]] = []
        gzip_members: list[dict[str, Any]] = []
        assignment_nnz = 0
        roles = self.manifest["archive_contract"]["roles"]

        for channel, accession in zip(
            self.channels,
            self.accessions,
            strict=True,
        ):
            modalities: dict[str, Any] = {}
            for modality, rows in self.rows.items():
                nnz = rows + channel
                matrix = {
                    "rows": rows,
                    "columns": self.columns,
                    "nnz": nnz,
                    "field": "integer",
                    "symmetry": "general",
                    "gzip_crc_verified": False,
                }
                features = {
                    "rows": rows,
                    "gzip_crc_verified": True,
                    "uncompressed_sha256": feature_hashes[modality],
                    "sample": [],
                }
                barcode = {
                    "rows": self.columns,
                    "gzip_crc_verified": True,
                    "uncompressed_sha256": barcode_hash,
                    "exact_axis_check": "uncompressed-byte identity",
                    "sample": [],
                }
                record: dict[str, Any] = {
                    "matrix": matrix,
                    "features": features,
                    "barcode": barcode,
                }
                if modality != "transcriptome":
                    relative = (
                        f"audit/support/channel{channel:02d}_"
                        f"{modality}_row_support.tsv.gz"
                    )
                    support_path = self.root / relative
                    support_counts = [1] * rows
                    support_counts[0] += nnz - rows
                    with gzip.open(
                        support_path,
                        mode="wt",
                        encoding="utf-8",
                        newline="",
                    ) as handle:
                        handle.write("row_index\tfeature_id\tfeature_name\tnnz_cells\n")
                        for index, count in enumerate(support_counts, start=1):
                            handle.write(
                                f"{index}\tfeature{index}\tfeature{index}\t{count}\n"
                            )
                    support_record = {
                        "path": relative,
                        "sha256": sha256(support_path),
                        "size_bytes": support_path.stat().st_size,
                    }
                    support_records.append(support_record)
                    record["assignment_scan"] = {
                        "active_cells": 3,
                        "coordinate_order": "column_major",
                        "declared_nnz": nnz,
                        "observed_nnz": nnz,
                        "gzip_crc_verified": True,
                        "supported_features": rows,
                        "row_support_sha256": hashlib.sha256(
                            array("Q", support_counts).tobytes()
                        ).hexdigest(),
                        "support_table": support_record,
                    }
                    assignment_nnz += nnz
                modalities[modality] = record

                prefix = f"{accession}_channel{channel}_{modality}_"
                for role, suffix in roles.items():
                    if role == "features":
                        member_hash = feature_hashes[modality]
                    elif role == "barcode":
                        member_hash = barcode_hash
                    else:
                        member_hash = hashlib.sha256(
                            f"{channel}-{modality}-matrix".encode()
                        ).hexdigest()
                    gzip_members.append(
                        {
                            "name": f"{prefix}{role}{suffix}",
                            "uncompressed_bytes": 10,
                            "uncompressed_sha256": member_hash,
                        }
                    )
            schema_channels.append(
                {
                    "channel": channel,
                    "accession": accession,
                    "barcode_axis_sha256": barcode_hash,
                    "barcode_axis": {
                        "rows": self.columns,
                        "expected_rows": self.columns,
                        "gzip_crc_verified": True,
                        "uncompressed_sha256": barcode_hash,
                        "uncompressed_byte_identity_across_modalities": True,
                        "sample": [],
                    },
                    "modalities": modalities,
                }
            )

        binding = self._binding()
        inventory = self._inventory()
        write_json(
            self.paths.support_report,
            {
                "status": "complete",
                "expected_support_tables": len(support_records),
                "support_tables": support_records,
                "input_binding": binding,
            },
        )
        write_json(
            self.paths.schema_report,
            {
                "status": "pass",
                "inventory": inventory,
                "feature_axis_sha256": feature_hashes,
                "barcode_axes_verified": True,
                "assignment_coordinates_scanned": True,
                "transcriptome_coordinates_scanned": False,
                "declared_assignment_nnz_total": assignment_nnz,
                "channels": schema_channels,
                "input_binding": binding,
                "support_completion_manifest": ("audit/support/COMPLETE.json"),
            },
        )
        write_json(
            self.paths.gzip_report,
            {
                "status": "pass",
                "inventory": inventory,
                "member_count": len(gzip_members),
                "all_members_read_to_gzip_eof": True,
                "members": gzip_members,
                "input_binding": binding,
            },
        )

    def _prepared_arrays(self) -> dict[str, np.ndarray]:
        totals = np.array([0, 50, 101, 120, 200, 500], dtype=np.int64)
        return {
            "totals": totals,
            "detected_features": np.array(
                [0, 2, 3, 3, 4, 4],
                dtype=np.uint16,
            ),
            "mitochondrial_counts": np.array(
                [0, 1, 2, 3, 4, 5],
                dtype=np.int64,
            ),
            "feature_totals": np.array(
                [300, 300, 200, 171],
                dtype=np.int64,
            ),
            "ambient_counts": np.array([20, 15, 10, 5], dtype=np.int64),
            "ambient_probability": np.array(
                [0.4, 0.3, 0.2, 0.1],
                dtype=np.float64,
            ),
            "tested_indices": np.array([2, 3, 4, 5], dtype=np.int32),
            "observed_data_log_probability": np.array(
                [-10.0, -11.0, -12.0, -13.0],
                dtype=np.float64,
            ),
        }

    def _write_channel(self, channel: int, accession: str) -> None:
        directory = self.paths.work / f"channel{channel:02d}"
        directory.mkdir()
        arrays_path = directory / "cell-call-inputs.npz"
        np.savez_compressed(arrays_path, **self._prepared_arrays())
        source_stage_sha = sha256(self.paths.source_stage_code)
        preprocess_sha = sha256(self.paths.preprocess_code)
        receipt_sha = sha256(self.paths.oracle_receipt)
        preflight_sha = sha256(self.paths.preflight_receipt)
        preflight_code_sha = sha256(self.paths.preflight_code)
        matrix_nnz = self.rows["transcriptome"] + channel
        prepare = {
            "status": "complete",
            "channel": channel,
            "accession": accession,
            "matrix_member": (
                f"{accession}_channel{channel}_transcriptome_matrix.mtx.gz"
            ),
            "raw_archive_sha256": sha256(self.paths.archive),
            "source_stage_sha256": source_stage_sha,
            "preprocess_sha256": preprocess_sha,
            "r_oracle_receipt_sha256": receipt_sha,
            "preflight_receipt_sha256": preflight_sha,
            "preflight_implementation_sha256": preflight_code_sha,
            "numpy_version": np.__version__,
            "scipy_version": source_artifacts.scipy.__version__,
            "arrays_path": (f"work/channel{channel:02d}/cell-call-inputs.npz"),
            "arrays_sha256": sha256(arrays_path),
            "shape": [self.rows["transcriptome"], self.columns],
            "nnz": matrix_nnz,
            "coordinate_order": "column_major",
            "nonzero_barcodes": 5,
            "tested_barcodes": 4,
            "maximum_total": 500,
            "knee": 150.0,
            "inflection": 100.0,
            "expressed_features": 4,
            "mitochondrial_features": 1,
            "ambient_barcodes": 2,
            "elapsed_seconds": 1.0,
        }
        prepare_path = directory / "prepare-report.json"
        write_json(prepare_path, prepare)

        tested = self._prepared_arrays()["tested_indices"]
        p_value = np.array(
            [1 / 10_001, 1_001 / 10_001, 5_001 / 10_001, 1.0],
            dtype=np.float64,
        )
        temporary = p_value.copy()
        temporary[np.array([False, False, True, True])] = 0.0
        fdr = bh(temporary)
        accepted = tested[fdr <= 0.01]
        limited = np.array([True, False, False, False], dtype=bool)
        output_path = directory / "cell-calls-production.npz"
        np.savez_compressed(
            output_path,
            tested_indices=tested,
            p_value=p_value,
            fdr=fdr,
            limited=limited,
            accepted_indices=accepted,
        )
        write_json(
            directory / "cell-calls-production.json",
            {
                "status": "complete",
                "stage": "production",
                "channel": channel,
                "iterations": 10_000,
                "iteration_block_size": 64,
                "prepare_report_sha256": sha256(prepare_path),
                "prepared_arrays_sha256": sha256(arrays_path),
                "source_stage_sha256": source_stage_sha,
                "preprocess_sha256": preprocess_sha,
                "r_oracle_receipt_sha256": receipt_sha,
                "preflight_receipt_sha256": preflight_sha,
                "preflight_implementation_sha256": preflight_code_sha,
                "numpy_version": np.__version__,
                "scipy_version": source_artifacts.scipy.__version__,
                "output_path": (f"work/channel{channel:02d}/cell-calls-production.npz"),
                "output_sha256": sha256(output_path),
                "tested_barcodes": tested.size,
                "accepted_barcodes": accepted.size,
                "forced_retained_barcodes": 2,
                "limited_barcodes": 1,
                "knee": 150.0,
                "elapsed_seconds": 1.0,
            },
        )

    def mutate_json(
        self,
        path: Path,
        mutation: Callable[[dict[str, Any]], None],
    ) -> None:
        payload = json.loads(path.read_text(encoding="utf-8"))
        mutation(payload)
        write_json(path, payload)

    def rewrite_npz(
        self,
        path: Path,
        mutation: Callable[[dict[str, np.ndarray]], None],
    ) -> None:
        with np.load(path, allow_pickle=False) as payload:
            arrays = {name: payload[name] for name in payload.files}
        mutation(arrays)
        np.savez_compressed(path, **arrays)


class SourceArtifactTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.fixture = SyntheticFixture(Path(self.temporary.name))

    def verify(self) -> dict[str, Any]:
        return source_artifacts.verify_source_artifacts(
            self.fixture.paths,
            self.fixture.contract,
        )

    def assert_fails(self, fragment: str) -> None:
        with self.assertRaisesRegex(
            source_artifacts.SourceArtifactError,
            fragment,
        ):
            self.verify()

    def test_passing_synthetic_fixture(self) -> None:
        result = self.verify()
        self.assertEqual(result["status"], "pass")
        self.assertEqual(result["total_tested_barcodes"], 8)
        self.assertEqual(result["total_accepted_barcodes"], 6)
        self.assertFalse(result["transcriptome_matrix_values_read"])

    def test_raw_archive_same_size_tamper_is_rejected(self) -> None:
        data = bytearray(self.fixture.paths.archive.read_bytes())
        data[0] ^= 1
        self.fixture.paths.archive.write_bytes(data)
        self.assert_fails("raw archive SHA-256 changed")

    def test_raw_archive_size_tamper_is_rejected(self) -> None:
        with self.fixture.paths.archive.open("ab") as handle:
            handle.write(b"x")
        self.assert_fails("raw archive byte count changed")

    def test_protocol_tamper_is_rejected_by_frozen_hash(self) -> None:
        self.fixture.mutate_json(
            self.fixture.paths.protocol,
            lambda payload: payload.update({"schema_version": 2}),
        )
        self.assert_fails("protocol SHA-256 changed")

    def test_manifest_tamper_is_rejected_by_frozen_hash(self) -> None:
        self.fixture.mutate_json(
            self.fixture.paths.manifest,
            lambda payload: payload.update({"dataset": "tampered"}),
        )
        self.assert_fails("manifest SHA-256 changed")

    def test_provenance_tamper_is_rejected_by_frozen_hash(self) -> None:
        self.fixture.mutate_json(
            self.fixture.paths.provenance,
            lambda payload: payload.update({"status": "tampered"}),
        )
        self.assert_fails("provenance SHA-256 changed")

    def test_each_audit_binding_tamper_is_rejected(self) -> None:
        for path in (
            self.fixture.paths.schema_report,
            self.fixture.paths.gzip_report,
            self.fixture.paths.support_report,
        ):
            with self.subTest(path=path.name):
                with tempfile.TemporaryDirectory() as directory:
                    fixture = SyntheticFixture(Path(directory))
                    target = fixture.paths.root / path.relative_to(
                        self.fixture.paths.root
                    )
                    fixture.mutate_json(
                        target,
                        lambda payload: payload["input_binding"].update(
                            {"manifest_sha256": "0" * 64}
                        ),
                    )
                    with self.assertRaisesRegex(
                        source_artifacts.SourceArtifactError,
                        "input binding mismatch",
                    ):
                        source_artifacts.verify_source_artifacts(
                            fixture.paths,
                            fixture.contract,
                        )

    def test_gzip_incomplete_report_is_rejected(self) -> None:
        self.fixture.mutate_json(
            self.fixture.paths.gzip_report,
            lambda payload: payload.update({"all_members_read_to_gzip_eof": False}),
        )
        self.assert_fails("did not read every member")

    def test_support_file_tamper_is_rejected(self) -> None:
        path = self.fixture.root / "audit/support/channel01_guides_row_support.tsv.gz"
        with path.open("ab") as handle:
            handle.write(b"tamper")
        self.assert_fails("byte count changed")

    def test_schema_support_binding_tamper_is_rejected(self) -> None:
        def mutate(payload: dict[str, Any]) -> None:
            scan = payload["channels"][0]["modalities"]["guides"]["assignment_scan"]
            scan["support_table"]["sha256"] = "0" * 64

        self.fixture.mutate_json(self.fixture.paths.schema_report, mutate)
        self.assert_fails("support bindings disagree")

    def test_coherently_rehashed_support_content_tamper_is_rejected(self) -> None:
        relative = "audit/support/channel01_guides_row_support.tsv.gz"
        path = self.fixture.root / relative
        with gzip.open(
            path,
            mode="wt",
            encoding="utf-8",
            newline="",
        ) as handle:
            handle.write("row_index\tfeature_id\tfeature_name\tnnz_cells\n")
            handle.write("1\tfeature1\tfeature1\t99\n")
            handle.write("2\tfeature2\tfeature2\t1\n")
            handle.write("3\tfeature3\tfeature3\t1\n")

        def update_support(payload: dict[str, Any]) -> None:
            record = next(
                item for item in payload["support_tables"] if item["path"] == relative
            )
            record["sha256"] = sha256(path)
            record["size_bytes"] = path.stat().st_size

        self.fixture.mutate_json(self.fixture.paths.support_report, update_support)

        def update_schema(payload: dict[str, Any]) -> None:
            record = payload["channels"][0]["modalities"]["guides"]["assignment_scan"][
                "support_table"
            ]
            record["sha256"] = sha256(path)
            record["size_bytes"] = path.stat().st_size

        self.fixture.mutate_json(self.fixture.paths.schema_report, update_schema)
        self.assert_fails("support nnz changed")

    def test_oracle_status_tamper_is_rejected(self) -> None:
        self.fixture.mutate_json(
            self.fixture.paths.oracle_receipt,
            lambda payload: payload.update({"status": "fail"}),
        )
        self.assert_fails("oracle receipt did not pass")

    def test_oracle_package_tamper_is_rejected(self) -> None:
        self.fixture.mutate_json(
            self.fixture.paths.oracle_receipt,
            lambda payload: payload.update({"DropletUtils_version": "0.0.0"}),
        )
        self.assert_fails("DropletUtils_version")

    def test_oracle_check_set_tamper_is_rejected(self) -> None:
        self.fixture.mutate_json(
            self.fixture.paths.oracle_receipt,
            lambda payload: payload["checks"].pop(),
        )
        self.assert_fails("check set changed")

    def test_source_preflight_status_tamper_is_rejected(self) -> None:
        self.fixture.mutate_json(
            self.fixture.paths.preflight_receipt,
            lambda payload: payload.update({"status": "fail"}),
        )
        self.assert_fails("source preflight receipt mismatch at status")

    def test_source_preflight_code_change_invalidates_receipt(self) -> None:
        self.fixture.paths.preflight_code.write_text(
            "def source_preflight():\n    return 108\n",
            encoding="utf-8",
        )
        self.assert_fails("preflight receipt mismatch at")

    def test_source_output_runtime_version_tamper_is_rejected(self) -> None:
        path = self.fixture.paths.work / "channel01/prepare-report.json"
        self.fixture.mutate_json(
            path,
            lambda payload: payload.update({"numpy_version": "0.0.0"}),
        )
        self.assert_fails("prepare report is stale at numpy_version")

    def test_preprocess_change_invalidates_oracle(self) -> None:
        self.fixture.paths.preprocess_code.write_text(
            "def preprocess():\n    return 108\n",
            encoding="utf-8",
        )
        self.assert_fails("source preflight implementation hashes changed")

    def test_r_oracle_change_invalidates_receipt(self) -> None:
        self.fixture.paths.r_oracle_code.write_text(
            "stopifnot(FALSE)\n",
            encoding="utf-8",
        )
        self.assert_fails("stale for r_oracle.R")

    def test_prepare_array_byte_tamper_is_rejected(self) -> None:
        path = self.fixture.paths.work / "channel01/cell-call-inputs.npz"
        data = bytearray(path.read_bytes())
        data[-1] ^= 1
        path.write_bytes(data)
        self.assert_fails("prepared arrays SHA-256 changed")

    def test_prepare_code_hash_tamper_is_rejected(self) -> None:
        path = self.fixture.paths.work / "channel01/prepare-report.json"
        self.fixture.mutate_json(
            path,
            lambda payload: payload.update({"source_stage_sha256": "0" * 64}),
        )
        self.assert_fails("prepare report is stale at source_stage_sha256")

    def test_source_stage_change_invalidates_prepare(self) -> None:
        self.fixture.paths.source_stage_code.write_text(
            "def source_stage():\n    return 108\n",
            encoding="utf-8",
        )
        self.assert_fails("source preflight implementation hashes changed")

    def test_reprepare_after_call_is_rejected(self) -> None:
        path = self.fixture.paths.work / "channel01/prepare-report.json"
        self.fixture.mutate_json(
            path,
            lambda payload: payload.update({"elapsed_seconds": 2.0}),
        )
        self.assert_fails("stale after prepare report replacement")

    def test_replaced_prepare_arrays_invalidate_call(self) -> None:
        directory = self.fixture.paths.work / "channel01"
        arrays_path = directory / "cell-call-inputs.npz"
        prepare_path = directory / "prepare-report.json"

        self.fixture.rewrite_npz(
            arrays_path,
            lambda arrays: arrays["observed_data_log_probability"].__setitem__(
                0,
                -10.5,
            ),
        )
        self.fixture.mutate_json(
            prepare_path,
            lambda payload: payload.update({"arrays_sha256": sha256(arrays_path)}),
        )
        self.assert_fails("stale after prepare report replacement")

    def test_call_report_wrong_iteration_count_is_rejected(self) -> None:
        path = self.fixture.paths.work / "channel01/cell-calls-production.json"
        self.fixture.mutate_json(
            path,
            lambda payload: payload.update({"iterations": 9_999}),
        )
        self.assert_fails("exactly 10000 production iterations")

    def test_call_output_byte_tamper_is_rejected(self) -> None:
        path = self.fixture.paths.work / "channel01/cell-calls-production.npz"
        data = bytearray(path.read_bytes())
        data[-1] ^= 1
        path.write_bytes(data)
        self.assert_fails("production call output SHA-256 changed")

    def test_call_tested_indices_tamper_is_rejected(self) -> None:
        directory = self.fixture.paths.work / "channel01"
        output_path = directory / "cell-calls-production.npz"
        report_path = directory / "cell-calls-production.json"
        self.fixture.rewrite_npz(
            output_path,
            lambda arrays: arrays.update(
                {
                    "tested_indices": np.array(
                        [2, 3, 4, 4],
                        dtype=np.int32,
                    )
                }
            ),
        )
        self.fixture.mutate_json(
            report_path,
            lambda payload: payload.update({"output_sha256": sha256(output_path)}),
        )
        self.assert_fails("not strictly increasing")

    def test_call_accepted_indices_tamper_is_rejected(self) -> None:
        directory = self.fixture.paths.work / "channel01"
        output_path = directory / "cell-calls-production.npz"
        report_path = directory / "cell-calls-production.json"

        def mutate(arrays: dict[str, np.ndarray]) -> None:
            arrays["accepted_indices"] = np.array([2, 3, 4, 5], dtype=np.int32)

        self.fixture.rewrite_npz(output_path, mutate)

        def update_report(payload: dict[str, Any]) -> None:
            payload["output_sha256"] = sha256(output_path)
            payload["accepted_barcodes"] = 4

        self.fixture.mutate_json(report_path, update_report)
        self.assert_fails("accepted_indices disagree")

    def test_call_fdr_tamper_is_rejected_even_when_rehashed(self) -> None:
        directory = self.fixture.paths.work / "channel01"
        output_path = directory / "cell-calls-production.npz"
        report_path = directory / "cell-calls-production.json"

        def mutate(arrays: dict[str, np.ndarray]) -> None:
            arrays["fdr"] = arrays["fdr"].copy()
            arrays["fdr"][1] = 0.123

        self.fixture.rewrite_npz(output_path, mutate)
        self.fixture.mutate_json(
            report_path,
            lambda payload: payload.update({"output_sha256": sha256(output_path)}),
        )
        self.assert_fails("do not reproduce frozen BH")

    def test_missing_channel_artifact_is_rejected(self) -> None:
        path = self.fixture.paths.work / "channel02/cell-calls-production.json"
        path.unlink()
        self.assert_fails("cannot resolve channel 2 production call report")

    def test_report_path_escape_is_rejected(self) -> None:
        path = self.fixture.paths.work / "channel01/prepare-report.json"
        self.fixture.mutate_json(
            path,
            lambda payload: payload.update({"arrays_path": "../../escape.npz"}),
        )
        self.assert_fails("prepare arrays path changed")

    def test_symlinked_support_table_is_rejected(self) -> None:
        path = self.fixture.root / "audit/support/channel01_guides_row_support.tsv.gz"
        target = self.fixture.root / "real-support.bin"
        target.write_bytes(path.read_bytes())
        path.unlink()
        path.symlink_to(target)
        self.assert_fails("cannot be a symlink")


if __name__ == "__main__":
    unittest.main()
