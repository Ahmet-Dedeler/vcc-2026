from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import source_stage  # noqa: E402


class SourceStageTests(unittest.TestCase):
    def test_atomic_npz_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "arrays.npz"
            source_stage.atomic_savez(
                path,
                a=np.arange(10, dtype=np.int64),
                b=np.array([True, False]),
            )
            with np.load(path, allow_pickle=False) as payload:
                np.testing.assert_array_equal(payload["a"], np.arange(10))
                np.testing.assert_array_equal(payload["b"], [True, False])

    def test_atomic_json_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report.json"
            source_stage.atomic_write_json(path, {"status": "pass", "value": 3})
            self.assertEqual(
                json.loads(path.read_text()),
                {"status": "pass", "value": 3},
            )

    def test_member_names_are_exact(self) -> None:
        record = {"accession": "GSM7951413", "channel": 1}
        self.assertEqual(
            source_stage.member_name(record, "transcriptome", "matrix"),
            "GSM7951413_channel1_transcriptome_matrix.mtx.gz",
        )
        self.assertEqual(
            source_stage.member_name(record, "guides", "features"),
            "GSM7951413_channel1_guides_features.tsv.gz",
        )


if __name__ == "__main__":
    unittest.main()
