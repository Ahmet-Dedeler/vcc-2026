#!/usr/bin/env python3
"""Small safety tests for verify_and_extract.py."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile


MODULE_PATH = Path(__file__).with_name("verify_and_extract.py")
SPEC = importlib.util.spec_from_file_location("verify_and_extract", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class VerifyAndExtractTests(unittest.TestCase):
    def test_accepts_regular_members(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive_path = root / "safe.zip"
            destination = root / "out"
            with ZipFile(archive_path, "w") as archive:
                archive.writestr("bundle/starter.toml", "[training]\n")

            with ZipFile(archive_path) as archive:
                members = MODULE.validate_members(archive, destination)

            self.assertEqual(members[0]["name"], "bundle/starter.toml")

    def test_rejects_parent_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive_path = root / "unsafe.zip"
            destination = root / "out"
            with ZipFile(archive_path, "w") as archive:
                archive.writestr("../escape.txt", "no")

            with ZipFile(archive_path) as archive:
                with self.assertRaisesRegex(ValueError, "unsafe ZIP member path"):
                    MODULE.validate_members(archive, destination)

    def test_unwraps_matching_single_top_level_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            staging = root / "staging"
            destination = root / "competition_support_set"
            wrapped = staging / destination.name
            wrapped.mkdir(parents=True)
            (wrapped / "starter.toml").write_text("[training]\n")

            self.assertEqual(
                MODULE.extracted_payload_root(staging, destination),
                wrapped,
            )

    def test_preserves_unwrapped_staging_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            staging = root / "staging"
            destination = root / "competition_support_set"
            staging.mkdir()
            (staging / "starter.toml").write_text("[training]\n")

            self.assertEqual(
                MODULE.extracted_payload_root(staging, destination),
                staging,
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
