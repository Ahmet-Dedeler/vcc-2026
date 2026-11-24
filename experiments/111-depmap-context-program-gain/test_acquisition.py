from __future__ import annotations

import copy
import csv
import hashlib
import io
import json
import stat
import struct
import sys
import tempfile
import unittest
import warnings
import zipfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import acquire  # noqa: E402


def make_zip(
    members: list[tuple[zipfile.ZipInfo | str, bytes]],
) -> bytes:
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for member, payload in members:
            archive.writestr(member, payload)
    return output.getvalue()


def archive_contract(
    payload: bytes,
    selected_names: list[str],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    tail_start = max(0, len(payload) - 131072)
    eocd = acquire.parse_eocd(payload[tail_start:], tail_start, len(payload))
    central = payload[eocd["offset"] : eocd["offset"] + eocd["size_bytes"]]
    entries = acquire.parse_central_directory(central, eocd["entry_count"])
    central_spec: dict[str, object] = {
        "entry_count": eocd["entry_count"],
        "offset": eocd["offset"],
        "size_bytes": eocd["size_bytes"],
        "sha256": hashlib.sha256(central).hexdigest(),
    }
    member_specs: list[dict[str, object]] = []
    for name in selected_names:
        entry = entries[name]
        member_specs.append(
            {
                "name": name,
                "compressed_size": entry.compressed_size,
                "uncompressed_size": entry.uncompressed_size,
                "crc32": f"{entry.crc32:08x}",
                "flag_bits": entry.flag_bits,
                "compression_method": entry.compression_method,
                "local_header_offset": entry.local_header_offset,
            }
        )
    return central_spec, member_specs


def write_csv(path: Path, rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle).writerows(rows)


class FakeResponse:
    def __init__(
        self,
        *,
        status: int,
        headers: dict[str, str],
        payload: bytes,
        url: str = "https://example.test/archive.zip",
    ) -> None:
        self.status = status
        self.headers = headers
        self.payload = payload
        self.url = url
        self.read_called = False
        self.closed = False

    def read(self, amount: int) -> bytes:
        self.read_called = True
        return self.payload[:amount]

    def geturl(self) -> str:
        return self.url

    def close(self) -> None:
        self.closed = True


class FakeHttpsClient(acquire.HttpsClient):
    def __init__(self, response: FakeResponse) -> None:
        self.response = response

    def _open(self, request: object) -> FakeResponse:
        del request
        return self.response


class ManifestContractTests(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = json.loads(acquire.DEFAULT_MANIFEST.read_text(encoding="utf-8"))
        self.protocol = json.loads(acquire.DEFAULT_PROTOCOL.read_text(encoding="utf-8"))

    def test_checked_in_contracts_and_zero_overlap_pass(self) -> None:
        report = acquire.run_manifest()
        self.assertEqual(report["cell_lines"], 16)
        self.assertEqual(report["lineages"], 11)
        self.assertEqual(report["targets"], 100)
        self.assertEqual(report["vcc_panel_overlap"], 0)
        self.assertFalse(report["selective_member_acquisition_authorized"])
        self.assertTrue(report["execution_amendment_selective_acquisition_authorized"])

    def test_cell_line_count_drift_fails(self) -> None:
        self.manifest["dataset_contract"]["cell_lines"].pop()
        with self.assertRaisesRegex(acquire.AcquisitionError, "exactly 16"):
            acquire.validate_manifest(self.manifest)

    def test_lineage_count_drift_fails(self) -> None:
        lines = self.manifest["dataset_contract"]["cell_lines"]
        old_lineage = lines[-1]["lineage"]
        self.assertEqual(
            sum(line["lineage"] == old_lineage for line in lines),
            1,
        )
        lines[-1]["lineage"] = lines[0]["lineage"]
        with self.assertRaisesRegex(acquire.AcquisitionError, "exactly 11"):
            acquire.validate_manifest(self.manifest)

    def test_target_count_drift_fails(self) -> None:
        self.manifest["dataset_contract"]["targets"].pop()
        with self.assertRaisesRegex(acquire.AcquisitionError, "100 unique targets"):
            acquire.validate_manifest(self.manifest)

    def test_nonzero_vcc_panel_overlap_fails(self) -> None:
        with acquire.REPO_ROOT.joinpath("data/challenge/pert_counts.csv").open(
            "r",
            encoding="utf-8",
            newline="",
        ) as handle:
            panel_target = next(csv.DictReader(handle))["target_gene"]
        self.manifest["dataset_contract"]["targets"][0] = panel_target
        acquire.validate_manifest(self.manifest)
        with self.assertRaisesRegex(acquire.AcquisitionError, "overlap drifted"):
            acquire.validate_challenge_panel(self.manifest)

    def test_exact_selective_member_allowlist_is_immutable(self) -> None:
        self.protocol["source"]["allowed_members"].append("processed/zscore_matrix.csv")
        with self.assertRaisesRegex(acquire.AcquisitionError, "allowlist drifted"):
            acquire.validate_protocol(self.protocol)

    def test_processed_archive_cannot_gain_a_full_sha256(self) -> None:
        processed = next(
            item
            for item in self.manifest["figshare"]["data_article"]["files"]
            if item["name"] == "processed.zip"
        )
        processed["sha256"] = "0" * 64
        with self.assertRaisesRegex(acquire.AcquisitionError, "SHA-256 status"):
            acquire.validate_manifest(self.manifest)

    def test_original_preregistration_remains_closed(self) -> None:
        self.assertFalse(
            self.protocol["current_authorization"]["selective_member_acquisition"]
        )

    def test_execution_authorization_tampering_fails(self) -> None:
        authorization = json.loads(
            acquire.DEFAULT_AUTHORIZATION.read_text(encoding="utf-8")
        )
        authorization["forbidden_actions"]["held_truth_access"] = True
        with self.assertRaisesRegex(
            acquire.AcquisitionError,
            "forbidden actions drifted",
        ):
            acquire.validate_execution_authorization(
                authorization,
                manifest_path=acquire.DEFAULT_MANIFEST,
                protocol_path=acquire.DEFAULT_PROTOCOL,
            )

    def test_dry_run_is_local_and_nonmutating(self) -> None:
        report = acquire.run_acquire(
            confirm_selective_download=False,
            dry_run=True,
        )
        self.assertEqual(report["selected_compressed_bytes"], 300011660)
        self.assertFalse(report["full_processed_archive_downloaded"])
        self.assertFalse(report["single_cell_data_downloaded"])
        self.assertTrue(report["execution_amendment_authorized"])
        self.assertFalse(report["original_preregistration_authorized"])


class HttpRangeTests(unittest.TestCase):
    def test_ignored_range_200_is_rejected_before_body_read(self) -> None:
        response = FakeResponse(
            status=200,
            headers={"Content-Length": "3233927690"},
            payload=b"must not be read",
        )
        client = FakeHttpsClient(response)
        with self.assertRaisesRegex(acquire.AcquisitionError, "expected 206"):
            client.read_range(
                "https://example.test/archive.zip",
                10,
                19,
                100,
            )
        self.assertFalse(response.read_called)
        self.assertTrue(response.closed)

    def test_content_range_total_must_match(self) -> None:
        response = FakeResponse(
            status=206,
            headers={
                "Content-Length": "10",
                "Content-Range": "bytes 10-19/101",
            },
            payload=b"0123456789",
        )
        client = FakeHttpsClient(response)
        with self.assertRaisesRegex(acquire.AcquisitionError, "Content-Range"):
            client.read_range(
                "https://example.test/archive.zip",
                10,
                19,
                100,
            )

    def test_content_encoding_is_rejected(self) -> None:
        response = FakeResponse(
            status=206,
            headers={
                "Content-Length": "10",
                "Content-Range": "bytes 10-19/100",
                "Content-Encoding": "gzip",
            },
            payload=b"0123456789",
        )
        client = FakeHttpsClient(response)
        with self.assertRaisesRegex(acquire.AcquisitionError, "Content-Encoding"):
            client.read_range(
                "https://example.test/archive.zip",
                10,
                19,
                100,
            )


class ZipContractTests(unittest.TestCase):
    def test_selective_deflate_extraction_reads_only_member_ranges(self) -> None:
        wanted = b"selected payload\n" * 200
        payload = make_zip(
            [
                ("processed/wanted.csv", wanted),
                ("processed/not-selected.csv", b"other payload" * 200),
            ]
        )
        central, specs = archive_contract(payload, ["processed/wanted.csv"])
        source = acquire.MemoryRangeSource(payload)
        _, selected = acquire.inspect_zip(
            source,
            central,
            specs,
            tail_probe_bytes=min(131072, len(payload)),
        )
        source.requests.clear()
        extracted, report = acquire.extract_member_bytes(
            source,
            selected["processed/wanted.csv"],
            max_chunk_bytes=17,
        )
        self.assertEqual(extracted, wanted)
        self.assertEqual(report["sha256"], hashlib.sha256(wanted).hexdigest())
        member = selected["processed/wanted.csv"]
        self.assertTrue(source.requests)
        self.assertTrue(
            all(
                member.data_offset <= start <= end <= member.data_end
                for start, end in source.requests
            )
        )

    def test_unsafe_path_is_rejected(self) -> None:
        payload = make_zip([("../escape.csv", b"x")])
        tail_start = max(0, len(payload) - 131072)
        eocd = acquire.parse_eocd(payload[tail_start:], tail_start, len(payload))
        central = payload[eocd["offset"] : eocd["offset"] + eocd["size_bytes"]]
        with self.assertRaisesRegex(acquire.AcquisitionError, "unsafe ZIP"):
            acquire.parse_central_directory(central, eocd["entry_count"])

    def test_duplicate_path_is_rejected(self) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            payload = make_zip([("same.csv", b"a"), ("same.csv", b"b")])
        tail_start = max(0, len(payload) - 131072)
        eocd = acquire.parse_eocd(payload[tail_start:], tail_start, len(payload))
        central = payload[eocd["offset"] : eocd["offset"] + eocd["size_bytes"]]
        with self.assertRaisesRegex(acquire.AcquisitionError, "duplicate ZIP"):
            acquire.parse_central_directory(central, eocd["entry_count"])

    def test_encrypted_flag_is_rejected(self) -> None:
        payload = bytearray(make_zip([("safe.csv", b"x")]))
        central, _ = archive_contract(bytes(payload), ["safe.csv"])
        flag_offset = int(central["offset"]) + 8
        current = struct.unpack_from("<H", payload, flag_offset)[0]
        struct.pack_into("<H", payload, flag_offset, current | 0x1)
        central_bytes = bytes(
            payload[
                int(central["offset"]) : int(central["offset"])
                + int(central["size_bytes"])
            ]
        )
        with self.assertRaisesRegex(acquire.AcquisitionError, "encrypted"):
            acquire.parse_central_directory(central_bytes, 1)

    def test_symlink_member_is_rejected(self) -> None:
        link = zipfile.ZipInfo("link")
        link.create_system = 3
        link.external_attr = (stat.S_IFLNK | 0o777) << 16
        payload = make_zip([(link, b"target")])
        tail_start = max(0, len(payload) - 131072)
        eocd = acquire.parse_eocd(payload[tail_start:], tail_start, len(payload))
        central = payload[eocd["offset"] : eocd["offset"] + eocd["size_bytes"]]
        with self.assertRaisesRegex(acquire.AcquisitionError, "symlink"):
            acquire.parse_central_directory(central, 1)

    def test_member_size_crc_and_offset_drift_fail(self) -> None:
        payload = make_zip([("processed/wanted.csv", b"payload")])
        central, specs = archive_contract(payload, ["processed/wanted.csv"])
        source = acquire.MemoryRangeSource(payload)
        tail_start = len(payload) - min(131072, len(payload))
        eocd = acquire.parse_eocd(payload[tail_start:], tail_start, len(payload))
        entries = acquire.parse_central_directory(
            payload[eocd["offset"] : eocd["offset"] + eocd["size_bytes"]],
            1,
        )
        del source, central
        for field, replacement in (
            ("compressed_size", int(specs[0]["compressed_size"]) + 1),
            ("uncompressed_size", int(specs[0]["uncompressed_size"]) + 1),
            ("crc32", "00000000"),
            ("local_header_offset", int(specs[0]["local_header_offset"]) + 1),
        ):
            with self.subTest(field=field):
                drifted = copy.deepcopy(specs)
                drifted[0][field] = replacement
                with self.assertRaisesRegex(
                    acquire.AcquisitionError,
                    "metadata drifted",
                ):
                    acquire.verify_selected_entries(entries, drifted)

    def test_corrupt_compressed_member_fails(self) -> None:
        payload = bytearray(make_zip([("processed/wanted.csv", b"A" * 1000)]))
        central, specs = archive_contract(bytes(payload), ["processed/wanted.csv"])
        clean_source = acquire.MemoryRangeSource(bytes(payload))
        _, selected = acquire.inspect_zip(
            clean_source,
            central,
            specs,
            tail_probe_bytes=min(131072, len(payload)),
        )
        member = selected["processed/wanted.csv"]
        payload[member.data_offset + 2] ^= 0xFF
        corrupt_source = acquire.MemoryRangeSource(bytes(payload))
        with self.assertRaises(acquire.AcquisitionError):
            acquire.extract_member_bytes(
                corrupt_source,
                member,
                max_chunk_bytes=32,
            )

    def test_normalized_code_tree_ignores_archive_wrapper_and_apple_junk(self) -> None:
        first = make_zip(
            [
                ("root-a/src/a.py", b"print('a')\n"),
                ("root-a/README.md", b"hello\n"),
            ]
        )
        second = make_zip(
            [
                ("root-b/README.md", b"hello\n"),
                ("__MACOSX/root-b/._README.md", b"junk"),
                ("root-b/.DS_Store", b"junk"),
                ("root-b/src/a.py", b"print('a')\n"),
            ]
        )
        self.assertEqual(
            acquire.normalized_code_tree_sha256(first)["sha256"],
            acquire.normalized_code_tree_sha256(second)["sha256"],
        )


class SchemaAndReceiptTests(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = acquire.validate_manifest(
            json.loads(acquire.DEFAULT_MANIFEST.read_text(encoding="utf-8"))
        )
        self.protocol = acquire.validate_protocol(
            json.loads(acquire.DEFAULT_PROTOCOL.read_text(encoding="utf-8"))
        )

    def test_cell_line_metadata_exact_contract(self) -> None:
        rows = [list(acquire.METADATA_COLUMNS)]
        rows.extend(
            [
                [
                    record["cell_line"],
                    record["arxspan_id"],
                    record["batch"],
                    record["lineage"],
                ]
                for record in self.manifest["dataset_contract"]["cell_lines"]
            ]
        )
        output = io.StringIO(newline="")
        csv.writer(output, lineterminator="\n").writerows(rows)
        report = acquire._validate_cell_line_metadata_bytes(
            output.getvalue().encode(),
            self.manifest,
        )
        self.assertEqual(report["row_count"], 16)
        self.assertEqual(report["lineage_count"], 11)

        rows.pop()
        output = io.StringIO(newline="")
        csv.writer(output, lineterminator="\n").writerows(rows)
        with self.assertRaisesRegex(acquire.AcquisitionError, "exactly 16"):
            acquire._validate_cell_line_metadata_bytes(
                output.getvalue().encode(),
                self.manifest,
            )

    def test_matrix_schema_rejects_axis_and_value_drift(self) -> None:
        expected_pairs = [("A", "T1"), ("A", "T2")]
        with tempfile.TemporaryDirectory(dir=HERE) as temporary:
            path = Path(temporary) / "matrix.csv"
            write_csv(
                path,
                [
                    ["", "A", "A"],
                    ["target", "T1", "T2"],
                    ["response_id", "", ""],
                    ["G1", "1", "2"],
                ],
            )
            report = acquire._validate_matrix_csv(
                path,
                expected_pairs,
                value_kind="count",
                minimum_response_rows=1,
            )
            self.assertEqual(report["row_count"], 1)

            write_csv(
                path,
                [
                    ["", "A", "A"],
                    ["target", "T2", "T1"],
                    ["G1", "1", "2"],
                ],
            )
            with self.assertRaisesRegex(acquire.AcquisitionError, "axis/order"):
                acquire._validate_matrix_csv(
                    path,
                    expected_pairs,
                    value_kind="count",
                    minimum_response_rows=1,
                )

            write_csv(
                path,
                [
                    ["", "A", "A"],
                    ["target", "T1", "T2"],
                    ["G1", "-1", "2"],
                ],
            )
            with self.assertRaisesRegex(acquire.AcquisitionError, "count"):
                acquire._validate_matrix_csv(
                    path,
                    expected_pairs,
                    value_kind="count",
                    minimum_response_rows=1,
                )

    def test_matrix_schema_accepts_pinned_empty_field_missingness(self) -> None:
        expected_pairs = [("A", "T1"), ("A", "T2")]
        with tempfile.TemporaryDirectory(dir=HERE) as temporary:
            path = Path(temporary) / "lfc.csv"
            write_csv(
                path,
                [
                    ["", "A", "A"],
                    ["grna_target", "T1", "T2"],
                    ["response_id", "", ""],
                    ["G1", "", "0.25"],
                ],
            )
            report = acquire._validate_matrix_csv(
                path,
                expected_pairs,
                value_kind="lfc",
                minimum_response_rows=1,
            )
            self.assertEqual(report["missing_value_count"], 1)

    def test_matrix_schema_rejects_textual_missingness(self) -> None:
        expected_pairs = [("A", "T1"), ("A", "T2")]
        with tempfile.TemporaryDirectory(dir=HERE) as temporary:
            path = Path(temporary) / "lfc.csv"
            write_csv(
                path,
                [
                    ["", "A", "A"],
                    ["grna_target", "T1", "T2"],
                    ["G1", "NA", "0.25"],
                ],
            )
            with self.assertRaisesRegex(acquire.AcquisitionError, "not numeric"):
                acquire._validate_matrix_csv(
                    path,
                    expected_pairs,
                    value_kind="lfc",
                    minimum_response_rows=1,
                )

    def test_extracted_tree_rejects_any_unlisted_file(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as temporary:
            root = Path(temporary) / "extracted"
            for name in acquire.ALLOWED_MEMBERS:
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"")
            (root / "processed/zscore_matrix.csv").write_bytes(b"")
            with self.assertRaisesRegex(acquire.AcquisitionError, "allowlist drifted"):
                acquire._audit_extracted_tree(root)

    def test_receipt_source_binding_tampering_fails_before_file_reuse(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as temporary:
            root = Path(temporary)
            receipt_path = root / "receipt.json"
            manifest_hash = acquire.sha256_file(acquire.DEFAULT_MANIFEST)
            protocol_hash = acquire.sha256_file(acquire.DEFAULT_PROTOCOL)
            binding = acquire.source_binding(
                self.manifest,
                self.protocol,
                manifest_hash,
                protocol_hash,
            )
            binding["figshare_data_article"]["files"][0]["md5"] = "0" * 32
            receipt_path.write_text(
                json.dumps(
                    {
                        "schema_version": acquire.RECEIPT_SCHEMA_VERSION,
                        "status": "complete",
                        "source_binding": binding,
                        "members": {},
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(acquire.AcquisitionError, "stale or tampered"):
                acquire.verify_receipt(
                    receipt_path,
                    root / "missing-extracted",
                    self.manifest,
                    self.protocol,
                    manifest_hash,
                    protocol_hash,
                    "0" * 64,
                    {
                        "status": "authorized_source_only_execution",
                        "authorized_actions": acquire.AUTHORIZED_ACTIONS,
                        "forbidden_actions": acquire.FORBIDDEN_ACTIONS,
                        "output_policy": {},
                    },
                )

    def test_small_archive_md5_and_sha256_are_both_required(self) -> None:
        payload = b"small archive"
        md5 = hashlib.md5(payload, usedforsecurity=False).hexdigest()
        sha256 = hashlib.sha256(payload).hexdigest()
        report = acquire._verify_small_archive(
            payload,
            len(payload),
            md5,
            sha256,
            "synthetic",
        )
        self.assertEqual(report["md5"], md5)
        with self.assertRaisesRegex(acquire.AcquisitionError, "MD5"):
            acquire._verify_small_archive(
                payload,
                len(payload),
                "0" * 32,
                sha256,
                "synthetic",
            )
        with self.assertRaisesRegex(acquire.AcquisitionError, "SHA-256"):
            acquire._verify_small_archive(
                payload,
                len(payload),
                md5,
                "0" * 64,
                "synthetic",
            )


if __name__ == "__main__":
    unittest.main()
