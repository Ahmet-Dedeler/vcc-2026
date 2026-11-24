from __future__ import annotations

import copy
import csv
import gzip
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import source_audit as audit  # noqa: E402


def gzip_csv(rows: list[list[str]]) -> bytes:
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", mtime=0) as compressed:
        with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
            writer = csv.writer(text, lineterminator="\n")
            writer.writerows(rows)
    return buffer.getvalue()


class ManifestTests(unittest.TestCase):
    def test_real_manifest_is_source_only_and_locked(self) -> None:
        manifest = audit.load_manifest()
        self.assertEqual(manifest["scope"], "source_only_metadata")
        self.assertEqual(
            set(manifest["authorization"]["eligible_large_assets"]),
            audit.LARGE_CANDIDATES,
        )
        self.assertEqual(manifest["authorization"]["authorized_large_assets"], [])
        self.assertFalse(manifest["authorization"]["download_geo_raw_archive"])
        self.assertFalse(manifest["assets"]["geo_raw_archive"]["download_authorized"])
        for name in audit.LARGE_CANDIDATES:
            self.assertFalse(audit.is_large_download_authorized(manifest, name))

    def test_only_low_lfc_assets_can_ever_be_authorized(self) -> None:
        manifest = audit.load_manifest()
        synthetic = copy.deepcopy(manifest)
        digest = "a" * 64
        synthetic["authorization"]["authorized_large_assets"] = ["low1_lfc_matrix"]
        synthetic["assets"]["low1_lfc_matrix"]["sha256"] = digest
        synthetic["assets"]["low1_lfc_matrix"]["download_authorized"] = True
        self.assertTrue(
            audit.is_large_download_authorized(synthetic, "low1_lfc_matrix")
        )
        self.assertFalse(
            audit.is_large_download_authorized(synthetic, "geo_raw_archive")
        )

    def test_missing_full_digest_fails_authorization(self) -> None:
        manifest = audit.load_manifest()
        synthetic = copy.deepcopy(manifest)
        synthetic["authorization"]["authorized_large_assets"] = ["low2_lfc_matrix"]
        synthetic["assets"]["low2_lfc_matrix"]["download_authorized"] = True
        self.assertFalse(
            audit.is_large_download_authorized(synthetic, "low2_lfc_matrix")
        )

    def test_unexpected_large_candidate_fails_closed(self) -> None:
        manifest = json.loads(audit.MANIFEST_PATH.read_text())
        manifest["authorization"]["eligible_large_assets"].append("raw")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(
                audit.SourceAuditError,
                "only Low1 and Low2",
            ):
                audit.load_manifest(path)


class ParserTests(unittest.TestCase):
    def test_bounded_gzip_header_parser(self) -> None:
        payload = gzip_csv(
            [
                ["gene", "ABC_P1", "DEF_P1", "DEF_P2"],
                ["G1", "0.1", "0.2", "0.3"],
            ]
        )
        header, labels = audit.parse_gzip_csv_header(payload)
        self.assertEqual(labels, ["ABC_P1", "DEF_P1", "DEF_P2"])
        self.assertTrue(header.endswith(b"\n"))

    def test_truncated_header_fails_closed(self) -> None:
        payload = gzip_csv([["gene", "ABC_P1"], ["G1", "0.1"]])
        with self.assertRaisesRegex(audit.SourceAuditError, "does not contain"):
            audit.parse_gzip_csv_header(payload[:12])

    def test_duplicate_header_label_fails_closed(self) -> None:
        payload = gzip_csv([["gene", "ABC_P1", "ABC_P1"]])
        with self.assertRaisesRegex(audit.SourceAuditError, "duplicated"):
            audit.parse_gzip_csv_header(payload)

    def test_panel_overlap_uses_exact_gene_prefix(self) -> None:
        panel = ["ABC", "AB", "DEF"]
        labels = ["ABC_P1", "ABC_P2", "ABCD_P1", "DEF_ENST1(W)"]
        overlap, multi = audit.panel_overlap(panel, labels)
        self.assertEqual(overlap, ["ABC", "DEF"])
        self.assertEqual(multi, {"ABC": ["ABC_P1", "ABC_P2"]})

    def test_guide_library_parser_handles_meng_controls(self) -> None:
        payload = gzip_csv(
            [
                ["AAAA", "GENE_P1_1", "crispr"],
                ["CCCC", "GENE_P1_2", "crispr"],
                ["GGGG", "NTC_501", "crispr"],
                ["TTTT", "NTC_MENG1", "crispr"],
            ]
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "library.csv.gz"
            path.write_bytes(payload)
            parsed = audit.parse_guide_library(path)
        self.assertEqual(parsed["row_count"], 4)
        self.assertEqual(
            parsed["promoter_labels"],
            ["GENE_P1", "NTC_501", "NTC_MENG1"],
        )

    def test_expected_value_drift_fails(self) -> None:
        with self.assertRaisesRegex(audit.SourceAuditError, "audit drift"):
            audit.require_expected({"count": 1}, {"count": 2})


class RealMetadataContractTests(unittest.TestCase):
    def test_frozen_small_metadata_and_offline_audit(self) -> None:
        manifest = audit.load_manifest()
        for asset in manifest["assets"].values():
            if asset["group"] == "metadata":
                audit.verify_small_asset(HERE / asset["local_path"], asset)
        report = audit.build_report(live_probe=False)
        self.assertEqual(
            report["status"],
            "source_metadata_verified_large_assets_locked",
        )
        self.assertEqual(report["summary"]["vcc_panel_gene_overlap"], 67)
        self.assertEqual(
            report["authorization"]["authorized_large_assets"],
            [],
        )
        self.assertFalse(report["interpretation"]["held_truth_opened"])
        self.assertFalse(report["interpretation"]["large_file_downloaded"])


if __name__ == "__main__":
    unittest.main()
