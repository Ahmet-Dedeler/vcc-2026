from __future__ import annotations

import gzip
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import acquire  # noqa: E402
import audit  # noqa: E402


def synthetic_xlsx(path: Path, rows: list[list[str]]) -> None:
    strings: list[str] = []
    indices: dict[str, int] = {}
    for row in rows:
        for value in row:
            if value not in indices:
                indices[value] = len(strings)
                strings.append(value)
    shared = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        f'<sst xmlns="{audit.XML_MAIN}" count="{sum(map(len, rows))}" '
        f'uniqueCount="{len(strings)}">'
        + "".join(f"<si><t>{value}</t></si>" for value in strings)
        + "</sst>"
    )
    sheet_rows = []
    for row_number, row in enumerate(rows, 1):
        cells = []
        for column_number, value in enumerate(row, 1):
            column = chr(ord("A") + column_number - 1)
            cells.append(
                f'<c r="{column}{row_number}" t="s"><v>{indices[value]}</v></c>'
            )
        sheet_rows.append(f'<row r="{row_number}">{"".join(cells)}</row>')
    sheet = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        f'<worksheet xmlns="{audit.XML_MAIN}"><sheetData>'
        + "".join(sheet_rows)
        + "</sheetData></worksheet>"
    )
    workbook = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        f'<workbook xmlns="{audit.XML_MAIN}" xmlns:r="{audit.XML_REL}">'
        '<sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets>'
        "</workbook>"
    )
    relationships = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        f'<Relationships xmlns="{audit.XML_PACKAGE_REL}">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        'Target="worksheets/sheet1.xml"/>'
        "</Relationships>"
    )
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", relationships)
        archive.writestr("xl/sharedStrings.xml", shared)
        archive.writestr("xl/worksheets/sheet1.xml", sheet)


class AcquisitionTests(unittest.TestCase):
    def test_geo_raw_download_remains_disabled_even_with_pseudobulk_token(self) -> None:
        manifest = acquire.load_manifest()
        with self.assertRaises(acquire.AcquisitionError):
            acquire.authorize_asset(
                manifest,
                "geo_raw_archive",
                "I_CONFIRM_GSE295214_PSEUDOBULK_DOWNLOAD",
            )

    def test_pseudobulk_download_without_token_fails(self) -> None:
        manifest = acquire.load_manifest()
        with self.assertRaisesRegex(acquire.AcquisitionError, "exact confirmation"):
            acquire.authorize_asset(manifest, "public_pseudobulk_database", None)

    def test_pseudobulk_download_is_authorized_by_exact_amendment(self) -> None:
        manifest = acquire.load_manifest()
        amendment = acquire.load_pseudobulk_amendment(manifest)
        asset = acquire.authorize_asset(
            manifest,
            "public_pseudobulk_database",
            amendment["authorization"]["confirmation_token"],
            amendment=amendment,
        )
        self.assertEqual(
            asset["sha256"], amendment["asset_binding"]["huggingface_lfs_oid_sha256"]
        )
        self.assertFalse(amendment["authorization"]["download_geo_raw_archive"])

    def test_metadata_is_authorized_without_confirmation(self) -> None:
        manifest = acquire.load_manifest()
        asset = acquire.authorize_asset(manifest, "supplementary_table_6", None)
        self.assertEqual(asset["group"], "metadata")

    def test_derived_metadata_cannot_be_downloaded(self) -> None:
        manifest = acquire.load_manifest()
        with self.assertRaisesRegex(acquire.AcquisitionError, "produced"):
            acquire.authorize_asset(manifest, "post_qc_target_axis", None)


class ParserTests(unittest.TestCase):
    def test_synthetic_library_workbook(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "library.xlsx"
            synthetic_xlsx(
                path,
                [
                    ["element", "sequence"],
                    ["GENE1_+_one_f_0_left", "AAAA"],
                    ["GENE1_+_one_f_1_left", "AAAA"],
                    ["GENE1_-_two_f_0_left", "CCCC"],
                    ["GENE2_+_one_f_0_left", "GGGG"],
                    ["non-targeting_00001_f_0_left", "TTTT"],
                ],
            )
            result = audit.parse_library_workbook(path)
        self.assertEqual(result["row_count"], 5)
        self.assertEqual(result["unique_oligo_sequences"], 4)
        self.assertEqual(result["unique_targeting_sequences"], 3)
        self.assertEqual(result["unique_negative_controls"], 1)
        self.assertEqual(result["literal_target_labels"], ["GENE1", "GENE2"])

    def test_unknown_library_element_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "library.xlsx"
            synthetic_xlsx(
                path,
                [["element", "sequence"], ["unparseable", "AAAA"]],
            )
            with self.assertRaisesRegex(acquire.AcquisitionError, "unrecognized"):
                audit.parse_library_workbook(path)

    def test_synthetic_geo_filelist(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "filelist.txt"
            path.write_text(
                "#Archive/File\tName\tTime\tSize\tType\n"
                "Archive\tGSE295214_RAW.tar\tnow\t300\tTAR\n"
                "File\tA.h5\tnow\t100\tH5\n"
                "File\tB.h5\tnow\t200\tH5\n"
            )
            result = audit.parse_geo_filelist(path)
        self.assertEqual(result["archive_size_bytes"], 300)
        self.assertEqual(result["h5_member_count"], 2)
        self.assertEqual(result["h5_member_bytes"], 300)

    def test_non_h5_archive_member_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "filelist.txt"
            path.write_text(
                "#Archive/File\tName\tTime\tSize\tType\n"
                "Archive\tGSE295214_RAW.tar\tnow\t10\tTAR\n"
                "File\tunexpected.csv\tnow\t10\tCSV\n"
            )
            with self.assertRaisesRegex(acquire.AcquisitionError, "non-H5"):
                audit.parse_geo_filelist(path)

    def test_synthetic_soft_inventory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "family.soft.gz"
            payload = "\n".join(
                [
                    "!Series_title = title",
                    "!Series_geo_accession = GSE295214",
                    "!Series_status = Public",
                    "!Series_submission_date = date",
                    "!Series_last_update_date = date",
                    "!Series_web_link = doi",
                    "!Series_overall_design = design",
                    "^SAMPLE = GSM1",
                    "!Sample_supplementary_file_1 = ftp://example/A.h5",
                    "^SAMPLE = GSM2",
                    "!Sample_supplementary_file_1 = ftp://example/B.h5",
                ]
            )
            with gzip.open(path, "wt") as handle:
                handle.write(payload)
            result = audit.parse_family_soft(path)
        self.assertEqual(result["sample_count"], "2")
        self.assertEqual(result["sample_h5_url_count"], "2")

    def test_duplicate_axis_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "axis.txt"
            path.write_text("A\nA\n")
            with self.assertRaisesRegex(acquire.AcquisitionError, "duplicates"):
                audit.read_axis(path)

    def test_expected_audit_drift_fails(self) -> None:
        with self.assertRaisesRegex(acquire.AcquisitionError, "audit drift"):
            audit.require_expected({"count": 2}, {"count": 3})


class RealMetadataContractTests(unittest.TestCase):
    def test_current_metadata_audit_matches_frozen_expectations(self) -> None:
        report = audit.build_report()
        self.assertEqual(
            report["status"],
            "metadata_audit_complete_large_download_locked",
        )
        self.assertFalse(report["interpretation"]["full_download_performed"])
        self.assertEqual(report["summary"]["vcc_panel_post_qc_overlap"], 64)
        self.assertEqual(
            report["summary"]["vcc_panel_library_targets_removed_by_qc"],
            ["PARP3", "SYVN1"],
        )

    def test_manifest_is_json_and_metadata_only(self) -> None:
        manifest = json.loads(acquire.MANIFEST_PATH.read_text())
        self.assertEqual(manifest["scope"], "metadata_only")
        self.assertEqual(
            manifest["official_sources"]["code"]["analysis_commit"],
            "aa848b9d517071d27b141e058e9bdae1966a869f",
        )
        self.assertEqual(
            manifest["official_sources"]["code"]["model_commit"],
            "29f44c98c5621d575a058b549e0fdde0ba31a730",
        )
        self.assertFalse(manifest["authorization"]["download_geo_raw_archive"])
        self.assertFalse(
            manifest["authorization"]["download_full_public_pseudobulk_database"]
        )
        amendment = acquire.load_pseudobulk_amendment(manifest)
        self.assertTrue(
            amendment["authorization"]["download_public_pseudobulk_database"]
        )
        self.assertFalse(amendment["authorization"]["download_geo_raw_archive"])


if __name__ == "__main__":
    unittest.main()
