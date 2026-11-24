from __future__ import annotations

import copy
import gzip
import io
import json
import sys
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import acquire_lfc_v1 as acquisition  # noqa: E402


class _RangeHandler(BaseHTTPRequestHandler):
    payload = b""
    last_modified = "Tue, 30 Jun 2026 21:16:28 GMT"
    ignore_range = False

    def log_message(self, format: str, *args: object) -> None:
        return

    def _headers(self, status: int, length: int, content_range: str | None) -> None:
        self.send_response(status)
        self.send_header("Content-Length", str(length))
        self.send_header("Last-Modified", self.last_modified)
        self.send_header("Content-Type", "application/x-gzip")
        self.send_header("Accept-Ranges", "bytes")
        if content_range is not None:
            self.send_header("Content-Range", content_range)
        self.end_headers()

    def do_HEAD(self) -> None:
        self._headers(200, len(self.payload), None)

    def do_GET(self) -> None:
        range_header = self.headers.get("Range")
        if range_header and not self.ignore_range:
            offset = int(range_header.removeprefix("bytes=").removesuffix("-"))
            body = self.payload[offset:]
            self._headers(
                206,
                len(body),
                f"bytes {offset}-{len(self.payload) - 1}/{len(self.payload)}",
            )
        else:
            body = self.payload
            self._headers(200, len(body), None)
        self.wfile.write(body)


def gzip_payload() -> bytes:
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", mtime=0) as handle:
        handle.write(b"gene,A_P1\nG1,0.25\n")
    return buffer.getvalue()


class AmendmentTests(unittest.TestCase):
    def test_real_amendment_authorizes_only_low1_low2(self) -> None:
        amendment = acquisition.load_amendment()
        self.assertEqual(tuple(amendment["assets"]), acquisition.ASSET_NAMES)
        self.assertFalse(
            amendment["authorization"]["download_geo_raw_archive"]
        )
        self.assertFalse(
            amendment["authorization"]["download_any_other_processed_asset"]
        )

    def test_added_asset_fails_closed(self) -> None:
        amendment = json.loads(acquisition.AMENDMENT_PATH.read_text())
        amendment["assets"]["raw"] = copy.deepcopy(
            amendment["assets"]["low1_lfc_matrix"]
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "amendment.json"
            path.write_text(json.dumps(amendment))
            with self.assertRaisesRegex(acquisition.AcquisitionError, "asset set"):
                acquisition.load_amendment(path)

    def test_destination_escape_fails_closed(self) -> None:
        amendment = json.loads(acquisition.AMENDMENT_PATH.read_text())
        amendment["assets"]["low1_lfc_matrix"]["destination"] = "/tmp/escape.gz"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "amendment.json"
            path.write_text(json.dumps(amendment))
            with self.assertRaisesRegex(acquisition.AcquisitionError, "escapes"):
                acquisition.load_amendment(path)


class TransferTests(unittest.TestCase):
    def setUp(self) -> None:
        _RangeHandler.payload = gzip_payload()
        _RangeHandler.ignore_range = False
        self.server = ThreadingHTTPServer(("127.0.0.1", 0), _RangeHandler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.url = f"http://127.0.0.1:{self.server.server_port}/asset.gz"

    def tearDown(self) -> None:
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()

    def asset(self) -> dict[str, object]:
        return {
            "filename": "asset.gz",
            "url": self.url,
            "size_bytes": len(_RangeHandler.payload),
            "last_modified_http": _RangeHandler.last_modified,
            "content_type": "application/x-gzip",
            "accept_ranges": "bytes",
            "etag": None,
        }

    def test_head_binding(self) -> None:
        result = acquisition.head_binding(self.asset())
        self.assertEqual(result["content_length"], len(_RangeHandler.payload))

    def test_resumable_download_uses_distinct_existing_prefix(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            partial = Path(directory) / ".asset.primary.part"
            partial.write_bytes(_RangeHandler.payload[:11])
            result = acquisition.download_resumable(
                self.asset(),
                partial,
                transfer_name="primary",
            )
            self.assertEqual(result["resumed_from_bytes"], 11)
            self.assertEqual(partial.read_bytes(), _RangeHandler.payload)

    def test_ignored_range_fails_closed_without_truncating_prefix(self) -> None:
        _RangeHandler.ignore_range = True
        with tempfile.TemporaryDirectory() as directory:
            partial = Path(directory) / ".asset.primary.part"
            prefix = _RangeHandler.payload[:11]
            partial.write_bytes(prefix)
            with self.assertRaisesRegex(acquisition.AcquisitionError, "HTTP 200"):
                acquisition.download_resumable(
                    self.asset(),
                    partial,
                    transfer_name="primary",
                )
            self.assertEqual(partial.read_bytes(), prefix)

    def test_gzip_and_bytewise_verification(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            left = Path(directory) / "left.gz"
            right = Path(directory) / "right.gz"
            left.write_bytes(_RangeHandler.payload)
            right.write_bytes(_RangeHandler.payload)
            self.assertTrue(acquisition.files_identical(left, right))
            report = acquisition.verify_complete_gzip(left)
            self.assertEqual(report["uncompressed_size_bytes"], 18)
            right.write_bytes(_RangeHandler.payload[:-1] + b"x")
            self.assertFalse(acquisition.files_identical(left, right))


if __name__ == "__main__":
    unittest.main()
