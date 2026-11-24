#!/usr/bin/env python3
"""Safety tests for prepare_support_view.py."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import h5py


MODULE_PATH = Path(__file__).with_name("prepare_support_view.py")
SPEC = importlib.util.spec_from_file_location("prepare_support_view", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PrepareSupportViewTests(unittest.TestCase):
    def test_adds_marker_without_changing_source(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.h5"
            clone = root / "clone.h5"
            with h5py.File(source, "w") as handle:
                handle.create_dataset("X", data=[1.0, 2.0])

            MODULE.clone_file(source, clone)
            self.assertTrue(MODULE.ensure_log1p_marker(clone))
            self.assertFalse(MODULE.ensure_log1p_marker(clone))

            with h5py.File(source, "r") as handle:
                self.assertNotIn("uns", handle)
            with h5py.File(clone, "r") as handle:
                marker = handle["uns/log1p"]
                self.assertEqual(marker.attrs["encoding-type"], "dict")
                self.assertEqual(marker.attrs["encoding-version"], "0.1.0")


if __name__ == "__main__":
    unittest.main(verbosity=2)
