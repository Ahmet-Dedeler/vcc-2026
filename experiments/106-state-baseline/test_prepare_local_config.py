#!/usr/bin/env python3
"""Safety tests for prepare_local_config.py."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


MODULE_PATH = Path(__file__).with_name("prepare_local_config.py")
SPEC = importlib.util.spec_from_file_location("prepare_local_config", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PrepareLocalConfigTests(unittest.TestCase):
    def test_rewrites_only_the_frozen_colab_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "starter.toml"
            support = root / "support"
            source.write_text(
                '[datasets]\n'
                'x = "/content/state/competition_support_set/{a,b}.h5"\n'
                '[training]\n'
                'x = "train"\n'
            )
            digest = __import__("hashlib").sha256(source.read_bytes()).hexdigest()

            with patch.object(MODULE, "EXPECTED_SOURCE_SHA256", digest):
                rendered = MODULE.render_local_config(source, support)

            self.assertIn(str(support.resolve()), rendered)
            self.assertNotIn(MODULE.COLAB_DATA_ROOT, rendered)
            self.assertIn('[training]\nx = "train"', rendered)

    def test_rejects_unexpected_source(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "starter.toml"
            source.write_text("[datasets]\n")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                MODULE.render_local_config(source, Path(tmp) / "support")


if __name__ == "__main__":
    unittest.main(verbosity=2)
