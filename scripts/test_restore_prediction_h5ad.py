"""Small, offline recovery regressions; no production artifacts are modified."""
import io
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest

import h5py
import numpy as np

from restore_prediction_h5ad import restore, storage_content_sha256, strings


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.source = self.root / "source.h5ad"
        with h5py.File(self.source, "w") as f:
            x = f.create_group("X")
            x.attrs["shape"] = [1002, 1]
            for name, values in (("data", np.ones(1002, dtype=np.float32)),
                                 ("indices", np.zeros(1002, dtype=np.int32)),
                                 ("indptr", np.arange(1003, dtype=np.int64))):
                x.create_dataset(name, data=values)
            obs = f.create_group("obs")
            for name, values in (("_index", [str(i) for i in range(1002)]),
                                 ("target_gene", ["T"] * 1002), ("context", ["A"] * 1002)):
                obs.create_dataset(name, data=values, dtype=h5py.string_dtype())
            f.create_group("var").create_dataset("_index", data=["G"], dtype=h5py.string_dtype())

    def tearDown(self):
        self.temp.cleanup()

    def package(self, *, bad=False, duplicate=False):
        path = self.root / "with ' quote; literal.vcc"
        data = b"broken zstd" if bad else subprocess.check_output(
            ["zstd", "-q", "-c", str(self.source)])
        with tarfile.open(path, "w") as archive:
            for _ in range(2 if duplicate else 1):
                member = tarfile.TarInfo("pred.h5ad.zst")
                member.size = len(data)
                archive.addfile(member, io.BytesIO(data))
        return path

    def test_semantic_round_trip_and_literal_paths(self):
        digest = storage_content_sha256(self.source)
        out = self.root / "restored ' file.h5ad"
        restore(self.package(), out, expected_content=digest)
        self.assertEqual(storage_content_sha256(out), digest)
        with h5py.File(out) as f:
            self.assertEqual(strings(f["obs"], "_index")[-1], "cell_1001")

    def test_raw_preserves_payload_bytes(self):
        out = self.root / "raw.h5ad"
        restore(self.package(), out, raw=True)
        self.assertEqual(out.read_bytes(), self.source.read_bytes())

    def test_wrong_hash_leaves_no_output(self):
        out = self.root / "out.h5ad"
        with self.assertRaisesRegex(ValueError, "content mismatch"):
            restore(self.package(), out, expected_content="0" * 64)
        self.assertFalse(out.exists())
        self.assertEqual(list(self.root.glob(".restore-*")), [])

    def test_existing_file_never_overwritten(self):
        with self.assertRaises(FileExistsError):
            restore(self.package(), self.source)

    def test_invalid_late_index_is_rejected(self):
        with h5py.File(self.source, "r+") as f:
            f["obs/_index"][1001] = "unexpected"
        with self.assertRaisesRegex(ValueError, "row 1001"):
            restore(self.package(), self.root / "out.h5ad")
        self.assertFalse((self.root / "out.h5ad").exists())

    def test_nullable_missing_index_rejected(self):
        with h5py.File(self.source, "r+") as f:
            del f["obs/_index"]
            group = f["obs"].create_group("_index")
            group.create_dataset("values", data=[str(i) for i in range(1002)],
                                 dtype=h5py.string_dtype())
            group.create_dataset("mask", data=[False] * 1001 + [True])
        with self.assertRaisesRegex(ValueError, "missing labels"):
            restore(self.package(), self.root / "out.h5ad")

    def test_broken_compression_rejected(self):
        with self.assertRaises((RuntimeError, BrokenPipeError)):
            restore(self.package(bad=True), self.root / "out.h5ad")
        self.assertFalse((self.root / "out.h5ad").exists())

    def test_duplicate_payload_rejected(self):
        with self.assertRaisesRegex(ValueError, "exactly one"):
            restore(self.package(duplicate=True), self.root / "out.h5ad")


if __name__ == "__main__":
    unittest.main()
