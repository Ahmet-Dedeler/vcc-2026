from __future__ import annotations

import gzip
import hashlib
import io
import json
import sys
import tarfile
import tempfile
import unittest
from array import array
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import acquire  # noqa: E402
import audit_archive as audit  # noqa: E402


def gzip_bytes(payload: bytes, *, mtime: int = 0) -> bytes:
    target = io.BytesIO()
    with gzip.GzipFile(fileobj=target, mode="wb", mtime=mtime) as stream:
        stream.write(payload)
    return target.getvalue()


def matrix_bytes(rows: int, columns: int, entries: list[tuple[int, int, int]]) -> bytes:
    lines = [
        "%%MatrixMarket matrix coordinate integer general",
        "% synthetic",
        f"{rows} {columns} {len(entries)}",
    ]
    lines.extend(f"{row} {column} {value}" for row, column, value in entries)
    return gzip_bytes(("\n".join(lines) + "\n").encode())


def matrix_bytes_with_declared_nnz(
    rows: int,
    columns: int,
    declared_nnz: int,
    entries: list[tuple[int, int, int]],
) -> bytes:
    lines = [
        "%%MatrixMarket matrix coordinate integer general",
        "% synthetic",
        f"{rows} {columns} {declared_nnz}",
    ]
    lines.extend(f"{row} {column} {value}" for row, column, value in entries)
    return gzip_bytes(("\n".join(lines) + "\n").encode())


class SyntheticDataset:
    def __init__(self, root: Path, *, channels: int = 2) -> None:
        self.root = root
        self.downloads = root / "downloads"
        self.downloads.mkdir()
        self.archive_path = self.downloads / "synthetic_RAW.tar"
        self.filelist_path = self.downloads / "filelist.txt"
        self.provenance_path = root / "PROVENANCE.json"
        self.manifest_path = root / "ACQUISITION_MANIFEST.json"
        self.guide_reference_path = self.downloads / "guides.fa.gz"
        self.supplementary_table_path = self.downloads / "supplement.xlsx"
        self.guide_reference_path.write_bytes(b">guide\nACGT\n")
        self.supplementary_table_path.write_bytes(b"synthetic workbook")
        self.members: dict[str, bytes] = {}
        channel_records = []
        modalities = {
            "transcriptome": {"feature_rows": 2, "barcode_columns": 3},
            "guides": {"feature_rows": 2, "barcode_columns": 3},
            "labels": {"feature_rows": 2, "barcode_columns": 3},
        }
        barcodes = gzip_bytes(b"BC1\nBC2\nBC3\n")
        for channel in range(1, channels + 1):
            accession = f"GSM{channel:07d}"
            channel_records.append({"channel": channel, "accession": accession})
            for modality in modalities:
                prefix = f"{accession}_channel{channel}_{modality}_"
                self.members[f"{prefix}barcode.tsv.gz"] = barcodes
                self.members[f"{prefix}features.tsv.gz"] = gzip_bytes(
                    f"{modality}_1\t{modality} one\n"
                    f"{modality}_2\t{modality} two\n".encode()
                )
                self.members[f"{prefix}matrix.mtx.gz"] = matrix_bytes(
                    2,
                    3,
                    [(1, 1, 2), (2, 2, 1), (1, 3, 1)],
                )
        self._write_archive()
        self._write_filelist()

        manifest = {
            "schema_version": 1,
            "dataset": "GSE249595",
            "assets": {
                "file_inventory": {
                    "group": "metadata",
                    "filename": "filelist.txt",
                    "url": "https://example.test/filelist.txt",
                    "size_bytes": self.filelist_path.stat().st_size,
                    "sha256": acquire.sha256_file(self.filelist_path),
                },
                "guide_reference": {
                    "group": "metadata",
                    "filename": "guides.fa.gz",
                    "url": "https://example.test/guides.fa.gz",
                    "size_bytes": self.guide_reference_path.stat().st_size,
                    "sha256": acquire.sha256_file(self.guide_reference_path),
                },
                "supplementary_table": {
                    "group": "metadata",
                    "filename": "supplement.xlsx",
                    "url": "https://example.test/supplement.xlsx",
                    "size_bytes": self.supplementary_table_path.stat().st_size,
                    "sha256": acquire.sha256_file(self.supplementary_table_path),
                },
                "raw_archive": {
                    "group": "raw",
                    "filename": self.archive_path.name,
                    "url": "https://example.test/synthetic_RAW.tar",
                    "size_bytes": self.archive_path.stat().st_size,
                    "sha256": None,
                    "sha256_status": "not_source_pinned",
                },
            },
            "archive_contract": {
                "expected_member_count": len(self.members),
                "channels": channel_records,
                "modalities": modalities,
                "roles": {
                    "barcode": ".tsv.gz",
                    "features": ".tsv.gz",
                    "matrix": ".mtx.gz",
                },
            },
            "safety": {
                "download_root": "downloads",
                "report_root": "audit",
                "max_assignment_nnz_total": 100,
                "assignment_coordinate_order": "row_or_column_major_strict",
                "never_scan_transcriptome_coordinates": True,
                "never_extract_archive": True,
            },
            "provenance_expectations": {
                "source.geo_series": "GSE249595",
                "screen.channels": channels,
                "status": "preregistered_before_acquisition",
            },
        }
        self.manifest = manifest
        self._write_manifest()
        self.provenance = {
            "source": {
                "geo_series": "GSE249595",
                "raw_archive_url": manifest["assets"]["raw_archive"]["url"],
                "guide_reference_url": manifest["assets"]["guide_reference"]["url"],
                "file_inventory_url": manifest["assets"]["file_inventory"]["url"],
                "supplementary_table_url": manifest["assets"]["supplementary_table"][
                    "url"
                ],
                "downloaded_metadata": {
                    asset["filename"]: {
                        "size_bytes": asset["size_bytes"],
                        "sha256": asset["sha256"],
                    }
                    for asset in manifest["assets"].values()
                    if asset["group"] == "metadata"
                },
            },
            "screen": {"channels": channels},
            "status": "preregistered_before_acquisition",
        }
        self._write_provenance()

    def _write_archive(self) -> None:
        with tarfile.open(self.archive_path, mode="w") as archive:
            for name, payload in sorted(self.members.items()):
                info = tarfile.TarInfo(name)
                info.size = len(payload)
                info.mtime = 0
                archive.addfile(info, io.BytesIO(payload))

    def _write_filelist(self) -> None:
        lines = [
            "#Archive/File\tName\tTime\tSize\tType",
            (
                f"Archive\t{self.archive_path.name}\t01/01/2026 00:00:00\t"
                f"{self.archive_path.stat().st_size}\tTAR"
            ),
        ]
        for name, payload in sorted(self.members.items()):
            file_type = "MTX" if name.endswith(".mtx.gz") else "TSV"
            lines.append(
                f"File\t{name}\t01/01/2026 00:00:00\t{len(payload)}\t{file_type}"
            )
        self.filelist_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def _write_manifest(self) -> None:
        self.manifest_path.write_text(json.dumps(self.manifest), encoding="utf-8")

    def _write_provenance(self) -> None:
        self.provenance_path.write_text(
            json.dumps(self.provenance),
            encoding="utf-8",
        )

    def refresh_archive_contract(self) -> None:
        self._write_archive()
        self._write_filelist()
        filelist = self.manifest["assets"]["file_inventory"]
        filelist["size_bytes"] = self.filelist_path.stat().st_size
        filelist["sha256"] = acquire.sha256_file(self.filelist_path)
        self.manifest["assets"]["raw_archive"]["size_bytes"] = (
            self.archive_path.stat().st_size
        )
        self.provenance["source"]["downloaded_metadata"]["filelist.txt"] = {
            "size_bytes": filelist["size_bytes"],
            "sha256": filelist["sha256"],
        }
        self.provenance["source"]["raw_archive_size_bytes"] = self.manifest["assets"][
            "raw_archive"
        ]["size_bytes"]
        self._write_manifest()
        self._write_provenance()

    def write_receipt(self, asset_names: set[str] | None = None) -> dict[str, object]:
        manifest_hash = acquire.sha256_file(self.manifest_path)
        names = asset_names or set(self.manifest["assets"])
        entries = {}
        for name in names:
            asset = self.manifest["assets"][name]
            local = acquire.verify_local_asset(
                self.downloads / asset["filename"],
                asset,
            )
            local["binding"] = acquire.asset_binding(name, asset, manifest_hash)
            entries[name] = local
        receipt: dict[str, object] = {
            "schema_version": acquire.RECEIPT_SCHEMA_VERSION,
            "dataset": self.manifest["dataset"],
            "manifest_sha256": manifest_hash,
            "assets": entries,
        }
        (self.downloads / acquire.RECEIPT_NAME).write_text(
            json.dumps(receipt),
            encoding="utf-8",
        )
        return receipt


class AcquisitionTests(unittest.TestCase):
    def test_raw_download_requires_explicit_confirmation(self) -> None:
        manifest = acquire.load_manifest()
        with self.assertRaises(acquire.AcquisitionError):
            acquire.select_assets(manifest, "raw", False)
        selected = acquire.select_assets(manifest, "raw", True)
        self.assertEqual([name for name, _ in selected], ["raw_archive"])

    def test_download_dry_run_writes_nothing(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as temporary:
            download_dir = Path(temporary) / "downloads"
            result = acquire.run_download(
                "metadata",
                download_dir=download_dir,
                dry_run=True,
            )
            self.assertTrue(result["dry_run"])
            self.assertFalse(download_dir.exists())
            self.assertTrue(
                all(entry["status"] == "dry_run" for entry in result["assets"].values())
            )

    def test_experiment_path_confinement(self) -> None:
        with self.assertRaises(acquire.AcquisitionError):
            acquire.require_under_experiment(Path("/tmp/outside-exp107"))

    def test_complete_partial_is_recovered_without_curl(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as temporary:
            download_dir = Path(temporary) / "downloads"
            download_dir.mkdir()
            payload = b"complete payload"
            asset = {
                "group": "metadata",
                "filename": "asset.bin",
                "url": "https://example.test/asset.bin",
                "size_bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
            manifest_hash = "a" * 64
            partial = download_dir / "asset.bin.part"
            partial.write_bytes(payload)
            acquire.atomic_write_json(
                acquire.resume_state_path(partial),
                {
                    "schema_version": acquire.RESUME_SCHEMA_VERSION,
                    "binding": acquire.asset_binding(
                        "asset",
                        asset,
                        manifest_hash,
                    ),
                    "remote": {
                        "etag": '"v1"',
                        "last_modified": None,
                        "content_length": len(payload),
                    },
                },
            )
            result = acquire.download_asset(
                "asset",
                asset,
                download_dir,
                manifest_sha256=manifest_hash,
                dry_run=False,
                probe=lambda _asset, _directory: {
                    "etag": '"v1"',
                    "last_modified": None,
                    "content_length": len(payload),
                },
            )
            self.assertEqual(result["status"], "recovered_source_verified_partial")
            self.assertFalse(partial.exists())
            self.assertFalse(acquire.resume_state_path(partial).exists())
            self.assertEqual((download_dir / "asset.bin").read_bytes(), payload)

    def test_final_and_partial_ambiguity_fails(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as temporary:
            download_dir = Path(temporary) / "downloads"
            download_dir.mkdir()
            payload = b"complete payload"
            asset = {
                "group": "metadata",
                "filename": "asset.bin",
                "url": "https://example.test/asset.bin",
                "size_bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
            (download_dir / "asset.bin").write_bytes(payload)
            (download_dir / "asset.bin.part").write_bytes(payload[:3])
            with self.assertRaisesRegex(
                acquire.AcquisitionError,
                "coexists with partial",
            ):
                acquire.download_asset(
                    "asset",
                    asset,
                    download_dir,
                    manifest_sha256="a" * 64,
                    dry_run=False,
                )

    def test_curl_command_is_transport_hardened(self) -> None:
        asset = {
            "url": "https://example.test/asset.bin",
            "size_bytes": 123,
        }
        command = acquire.curl_command(
            asset,
            Path("/tmp/asset.part"),
            Path("/tmp/headers"),
            if_range='"v1"',
        )
        self.assertEqual(command[:2], ["curl", "--disable"])
        self.assertIn("--proto", command)
        self.assertIn("--proto-redir", command)
        self.assertIn("--max-filesize", command)
        self.assertIn('If-Range: "v1"', command)

    def test_manifest_cannot_reclassify_raw_as_metadata(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as temporary:
            path = Path(temporary) / "manifest.json"
            manifest = json.loads(acquire.DEFAULT_MANIFEST.read_text(encoding="utf-8"))
            manifest["assets"]["raw_archive"]["group"] = "metadata"
            manifest["assets"]["raw_archive"]["sha256"] = "0" * 64
            path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(
                acquire.AcquisitionError,
                "asset group must be",
            ):
                acquire.load_manifest(path)

    def test_manifest_member_count_mismatch_fails(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as temporary:
            path = Path(temporary) / "manifest.json"
            manifest = json.loads(acquire.DEFAULT_MANIFEST.read_text(encoding="utf-8"))
            manifest["archive_contract"]["expected_member_count"] = 143
            path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(
                acquire.AcquisitionError,
                "expected_member_count",
            ):
                acquire.load_manifest(path)

    def test_manifest_requires_all_safety_guards(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as temporary:
            path = Path(temporary) / "manifest.json"
            manifest = json.loads(acquire.DEFAULT_MANIFEST.read_text(encoding="utf-8"))
            del manifest["safety"]["never_extract_archive"]
            path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(
                acquire.AcquisitionError,
                "never_extract_archive",
            ):
                acquire.load_manifest(path)

    def test_changed_manifest_drops_stale_raw_receipt(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as temporary:
            dataset = SyntheticDataset(Path(temporary))
            dataset.write_receipt()
            dataset.manifest["assets"]["raw_archive"]["url"] = (
                "https://example.test/changed_RAW.tar"
            )
            dataset.provenance["source"]["raw_archive_url"] = dataset.manifest[
                "assets"
            ]["raw_archive"]["url"]
            dataset._write_manifest()
            dataset._write_provenance()

            acquire.run_download(
                "metadata",
                manifest_path=dataset.manifest_path,
                download_dir=dataset.downloads,
            )
            receipt = acquire.load_receipt(dataset.downloads)
            self.assertNotIn("raw_archive", receipt["assets"])
            with self.assertRaisesRegex(
                acquire.AcquisitionError,
                "has no receipt bound",
            ):
                acquire.download_asset(
                    "raw_archive",
                    dataset.manifest["assets"]["raw_archive"],
                    dataset.downloads,
                    manifest_sha256=acquire.sha256_file(dataset.manifest_path),
                    dry_run=False,
                )

    def test_raw_resume_requires_matching_remote_validator(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as temporary:
            dataset = SyntheticDataset(Path(temporary))
            asset = dataset.manifest["assets"]["raw_archive"]
            partial = dataset.downloads / f"{asset['filename']}.part"
            dataset.archive_path.replace(partial)
            acquire.atomic_write_json(
                acquire.resume_state_path(partial),
                {
                    "schema_version": acquire.RESUME_SCHEMA_VERSION,
                    "binding": acquire.asset_binding(
                        "raw_archive",
                        asset,
                        acquire.sha256_file(dataset.manifest_path),
                    ),
                    "remote": {
                        "etag": '"old"',
                        "last_modified": None,
                        "content_length": asset["size_bytes"],
                    },
                },
            )
            with self.assertRaisesRegex(
                acquire.AcquisitionError,
                "remote validators changed",
            ):
                acquire.download_asset(
                    "raw_archive",
                    asset,
                    dataset.downloads,
                    manifest_sha256=acquire.sha256_file(dataset.manifest_path),
                    dry_run=False,
                    probe=lambda _asset, _directory: {
                        "etag": '"new"',
                        "last_modified": None,
                        "content_length": asset["size_bytes"],
                    },
                )


class ArchiveAuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(dir=HERE)
        self.root = Path(self.temporary.name)
        self.dataset = SyntheticDataset(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_inventory_and_sparse_schema_scan(self) -> None:
        inventory, members = audit.inventory_archive(
            self.dataset.archive_path,
            self.dataset.filelist_path,
            self.dataset.manifest,
        )
        self.assertEqual(inventory["member_count"], 18)
        self.assertEqual(len(members), 18)

        support_dir = self.root / "support"
        report = audit.schema_audit(
            self.dataset.archive_path,
            self.dataset.filelist_path,
            self.dataset.manifest,
            verify_barcode_hashes=True,
            scan_assignments=True,
            support_dir=support_dir,
            input_binding={"raw_archive_sha256": "synthetic"},
        )
        self.assertTrue(report["barcode_axes_verified"])
        self.assertTrue(report["assignment_coordinates_scanned"])
        self.assertFalse(report["transcriptome_coordinates_scanned"])
        self.assertEqual(report["declared_assignment_nnz_total"], 12)
        guide_scan = report["channels"][0]["modalities"]["guides"]["assignment_scan"]
        self.assertEqual(guide_scan["observed_nnz"], 3)
        self.assertEqual(guide_scan["supported_features"], 2)
        self.assertEqual(guide_scan["active_cells"], 3)
        self.assertTrue((support_dir / "channel01_guides_row_support.tsv.gz").is_file())
        self.assertTrue((support_dir / "COMPLETE.json").is_file())
        completion = json.loads((support_dir / "COMPLETE.json").read_text())
        self.assertEqual(
            completion["input_binding"]["raw_archive_sha256"],
            "synthetic",
        )

    def test_missing_member_fails_closed(self) -> None:
        missing_name = next(iter(self.dataset.members))
        del self.dataset.members[missing_name]
        self.dataset._write_archive()
        with self.assertRaisesRegex(audit.AuditError, "archive member mismatch"):
            audit.inventory_archive(
                self.dataset.archive_path,
                self.dataset.filelist_path,
                self.dataset.manifest,
            )

    def test_matrix_dimension_mismatch_fails(self) -> None:
        matrix_name = next(
            name for name in self.dataset.members if "_guides_matrix.mtx.gz" in name
        )
        self.dataset.members[matrix_name] = matrix_bytes(
            3,
            3,
            [(1, 1, 1)],
        )
        self.dataset._write_archive()
        self.dataset._write_filelist()
        self.dataset.manifest["assets"]["raw_archive"]["size_bytes"] = (
            self.dataset.archive_path.stat().st_size
        )
        with self.assertRaisesRegex(audit.AuditError, "feature dimension"):
            audit.schema_audit(
                self.dataset.archive_path,
                self.dataset.filelist_path,
                self.dataset.manifest,
            )

    def test_barcode_row_count_mismatch_fails(self) -> None:
        for name in list(self.dataset.members):
            if name.endswith("_barcode.tsv.gz"):
                self.dataset.members[name] = gzip_bytes(b"ONLY_ONE\n")
        self.dataset.refresh_archive_contract()
        with self.assertRaisesRegex(audit.AuditError, "barcode row count"):
            audit.schema_audit(
                self.dataset.archive_path,
                self.dataset.filelist_path,
                self.dataset.manifest,
                verify_barcode_hashes=True,
            )

    def test_barcode_identity_uses_uncompressed_bytes(self) -> None:
        barcode_name = next(
            name for name in self.dataset.members if "_guides_barcode.tsv.gz" in name
        )
        self.dataset.members[barcode_name] = gzip_bytes(
            b"BC1\nBC2\nBC3\n",
            mtime=1,
        )
        self.dataset.refresh_archive_contract()
        report = audit.schema_audit(
            self.dataset.archive_path,
            self.dataset.filelist_path,
            self.dataset.manifest,
            verify_barcode_hashes=True,
        )
        self.assertTrue(report["barcode_axes_verified"])
        channel = report["channels"][0]
        self.assertTrue(
            channel["barcode_axis"]["uncompressed_byte_identity_across_modalities"]
        )

    def test_barcode_content_mismatch_fails(self) -> None:
        barcode_name = next(
            name for name in self.dataset.members if "_guides_barcode.tsv.gz" in name
        )
        self.dataset.members[barcode_name] = gzip_bytes(b"BC1\nBCX\nBC3\n")
        self.dataset.refresh_archive_contract()
        with self.assertRaisesRegex(
            audit.AuditError,
            "barcode files are not byte-identical",
        ):
            audit.schema_audit(
                self.dataset.archive_path,
                self.dataset.filelist_path,
                self.dataset.manifest,
                verify_barcode_hashes=True,
            )

    def test_extra_assignment_entry_aborts_at_first_excess_line(self) -> None:
        matrix_name = next(
            name for name in self.dataset.members if "_guides_matrix.mtx.gz" in name
        )
        self.dataset.members[matrix_name] = matrix_bytes_with_declared_nnz(
            2,
            3,
            1,
            [(1, 1, 1), (2, 2, 1), (1, 3, 1)],
        )
        self.dataset.refresh_archive_contract()
        with self.assertRaisesRegex(
            audit.AuditError,
            "more entries than declared at line 2",
        ):
            audit.schema_audit(
                self.dataset.archive_path,
                self.dataset.filelist_path,
                self.dataset.manifest,
                scan_assignments=True,
                support_dir=self.root / "support",
            )

    def test_negative_assignment_value_fails(self) -> None:
        matrix_name = next(
            name for name in self.dataset.members if "_guides_matrix.mtx.gz" in name
        )
        self.dataset.members[matrix_name] = matrix_bytes(2, 3, [(1, 1, -5)])
        self.dataset.refresh_archive_contract()
        with self.assertRaisesRegex(audit.AuditError, "must be positive"):
            audit.schema_audit(
                self.dataset.archive_path,
                self.dataset.filelist_path,
                self.dataset.manifest,
                scan_assignments=True,
                support_dir=self.root / "support",
            )

    def test_duplicate_assignment_coordinate_fails(self) -> None:
        matrix_name = next(
            name for name in self.dataset.members if "_guides_matrix.mtx.gz" in name
        )
        self.dataset.members[matrix_name] = matrix_bytes(
            2,
            3,
            [(1, 1, 1), (1, 1, 2)],
        )
        self.dataset.refresh_archive_contract()
        with self.assertRaisesRegex(audit.AuditError, "duplicate coordinate"):
            audit.schema_audit(
                self.dataset.archive_path,
                self.dataset.filelist_path,
                self.dataset.manifest,
                scan_assignments=True,
                support_dir=self.root / "support",
            )

    def test_failed_scan_cleans_staged_outputs(self) -> None:
        matrix_names = sorted(
            name for name in self.dataset.members if "_guides_matrix.mtx.gz" in name
        )
        self.dataset.members[matrix_names[-1]] = matrix_bytes(
            2,
            3,
            [(1, 1, -1)],
        )
        self.dataset.refresh_archive_contract()
        support_dir = self.root / "support"
        support_dir.mkdir()
        sentinel = support_dir / "previous-complete.txt"
        sentinel.write_text("preserve me", encoding="utf-8")
        with self.assertRaisesRegex(audit.AuditError, "must be positive"):
            audit.schema_audit(
                self.dataset.archive_path,
                self.dataset.filelist_path,
                self.dataset.manifest,
                scan_assignments=True,
                support_dir=support_dir,
            )
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "preserve me")
        self.assertEqual(list(support_dir.glob("*_row_support.tsv.gz")), [])
        self.assertEqual(
            list(self.root.glob(".support.stage-*")),
            [],
        )

    def test_path_traversal_member_fails(self) -> None:
        with tarfile.open(self.dataset.archive_path, mode="a") as archive:
            info = tarfile.TarInfo("../escape.tsv.gz")
            payload = gzip_bytes(b"x\n")
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
        self.dataset._write_filelist()
        self.dataset.manifest["assets"]["raw_archive"]["size_bytes"] = (
            self.dataset.archive_path.stat().st_size
        )
        with self.assertRaisesRegex(audit.AuditError, "unsafe"):
            audit.inventory_archive(
                self.dataset.archive_path,
                self.dataset.filelist_path,
                self.dataset.manifest,
            )

    def test_provenance_mismatch_fails(self) -> None:
        provenance = json.loads(
            self.dataset.provenance_path.read_text(encoding="utf-8")
        )
        provenance["screen"]["channels"] = 99
        self.dataset.provenance_path.write_text(
            json.dumps(provenance),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(audit.AuditError, "provenance mismatch"):
            audit.verify_provenance(
                self.dataset.manifest_path,
                self.dataset.provenance_path,
                self.dataset.downloads,
            )

    def test_same_size_raw_tampering_fails_required_provenance(self) -> None:
        self.dataset.write_receipt({"raw_archive"})
        with self.dataset.archive_path.open("r+b") as handle:
            handle.seek(17)
            original = handle.read(1)
            handle.seek(17)
            handle.write(bytes([original[0] ^ 1]))
        with self.assertRaisesRegex(audit.AuditError, "receipt mismatch"):
            audit.verify_provenance(
                self.dataset.manifest_path,
                self.dataset.provenance_path,
                self.dataset.downloads,
                require_raw=True,
            )

    def test_verify_raw_hash_requires_raw(self) -> None:
        with self.assertRaisesRegex(audit.AuditError, "requires --require-raw"):
            audit.verify_provenance(
                self.dataset.manifest_path,
                self.dataset.provenance_path,
                self.dataset.downloads,
                verify_raw_hash=True,
            )

    def test_stale_report_is_invalidated_before_audit(self) -> None:
        report_path = self.root / "stale-report.json"
        report_path.write_text('{"status":"pass"}', encoding="utf-8")
        resolved = audit.invalidate_report(report_path)
        self.assertEqual(resolved, report_path.resolve())
        self.assertFalse(report_path.exists())

    def test_gzip_integrity_reads_every_member(self) -> None:
        report = audit.verify_all_gzip_members(
            self.dataset.archive_path,
            self.dataset.filelist_path,
            self.dataset.manifest,
        )
        self.assertEqual(report["member_count"], 18)
        self.assertTrue(report["all_members_read_to_gzip_eof"])

    def test_gzip_integrity_rejects_bad_crc(self) -> None:
        member_name = next(iter(self.dataset.members))
        payload = bytearray(self.dataset.members[member_name])
        payload[-8] ^= 1
        self.dataset.members[member_name] = bytes(payload)
        self.dataset.refresh_archive_contract()
        with self.assertRaisesRegex(audit.AuditError, "gzip integrity audit failed"):
            audit.verify_all_gzip_members(
                self.dataset.archive_path,
                self.dataset.filelist_path,
                self.dataset.manifest,
            )

    def test_support_hash_is_stable(self) -> None:
        counts = array("Q", [2, 1])
        first = self.root / "first.tsv.gz"
        second = self.root / "second.tsv.gz"
        rows = [["g1", "Gene 1"], ["g2", "Gene 2"]]
        first_report = audit.write_support_table(first, rows, counts)
        second_report = audit.write_support_table(second, rows, counts)
        with gzip.open(first, "rb") as handle:
            first_payload = handle.read()
        with gzip.open(second, "rb") as handle:
            second_payload = handle.read()
        self.assertEqual(
            hashlib.sha256(first_payload).hexdigest(),
            hashlib.sha256(second_payload).hexdigest(),
        )
        self.assertEqual(first_report["sha256"], second_report["sha256"])


if __name__ == "__main__":
    unittest.main()
