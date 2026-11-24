"""Independent guard regression tests: all expression/releases/markers are synthetic.

No production source is opened, even for a digest. Fixtures live in a temporary
directory and patch every path; h5py value reads are instrumented. No repository
output is written. A separate real-input schema smoke test reads ONLY the saved
symmetry302.csv, COHORT.tsv and R5-ROWS.tsv text; HDF5/access writes are blocked.
Run:

    PYTHONDONTWRITEBYTECODE=1 uv run python test_access_guard303.py

Coordinator receipt mode (stdout only, never creates a release/marker):

    ... test_access_guard303.py --manifest /path/to/frozen-core-manifest.json

This binds this guard suite to a core manifest; it does not assert that analysis
tests ran. A separate frozen_analysis303 receipt remains mandatory. The code
hashes in that manifest must match the executing access/lib/test sources.
"""
from __future__ import annotations

import argparse
import copy
import csv
import io
import json
import os
import shutil
import stat
import sys
import tempfile
import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import access303 as A  # noqa: E402

L = A.L
REAL_HERE = Path(__file__).resolve().parent
REAL_SYM302 = L.SYM302
REAL_WORK = L.WORK


class AccessGuardTests(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.root = Path(self.stack.enter_context(
            tempfile.TemporaryDirectory(prefix="synthetic-access303-"))).resolve()
        real_root = L.ROOT
        for name in ("ROOT", "HERE", "WORK", "DB", "H1BANK", "CANON", "HELD279",
                     "ZNF279", "PERSRC290", "SYM302", "SUMS302"):
            old = getattr(L, name)
            self.stack.enter_context(patch.object(L, name, self.root / old.relative_to(real_root)))
        for name in ("RELEASE", "MARKER", "PRODUCT", "PREREG_FINAL"):
            old = getattr(A, name)
            self.stack.enter_context(patch.object(A, name, self.root / old.relative_to(real_root)))
        self.stack.enter_context(patch.object(L, "BLOCK", 2))
        L.WORK.mkdir(parents=True)
        L.H1BANK.mkdir(parents=True)
        L.DB.parent.mkdir(parents=True)
        L.CANON.parent.mkdir(parents=True)
        for name in A.CODE_FILES:
            shutil.copyfile(REAL_HERE / name, L.HERE / name)
        self.genes = ["A", "B", "C", "D"]
        self.rows = np.asarray([1, 2, 3], np.int64)
        self.cpm = np.asarray([[2500, 2500, 2500, 2500],
                               [3000, 2000, 2500, 2500],
                               [2000, 3000, 2500, 2500],
                               [2500, 2500, 3000, 2000]], np.float64)
        self.control = np.full((4, 4), np.log1p(2500), np.float64)
        with h5py.File(L.DB, "w") as f:
            f["X"] = np.log1p(self.cpm) - self.control
            f["layers/control"] = self.control
            f["layers/measure_mask"] = np.ones((4, 4), np.float32)

            def cat(name, values):
                cats, codes = np.unique(values, return_inverse=True)
                f[name + "/categories"] = np.asarray(cats, dtype="S")
                f[name + "/codes"] = codes.astype(np.int32)

            cat("obs/dataset", [L.BANK_DATASET] + [L.RESERVE_DATASET] * 3)
            cat("obs/perturbed_gene_name", ["Z", "A", "B", "C"])
            cat("obs/perturbed_gene", ["ENSGZ", "ENSGA", "ENSGB", "ENSGC"])
            cat("obs/condition", ["KD"] + ["hESC_TF_screen"] * 3)
            f["obs/id"] = np.asarray(["bank_Z", "hESC_TF_screen_A",
                                      "hESC_TF_screen_B", "hESC_TF_screen_C"], dtype="S")
            f["obs/perturbation_sign"] = np.asarray([-1, 1, 1, 1], np.int64)
            f["var/human_id"] = np.asarray(["ENSG" + x for x in self.genes], dtype="S")
            cat("var/gene_name", self.genes)
        L.CANON.write_text("gene_name\n" + "\n".join(self.genes) + "\n")
        (L.H1BANK / "genes.txt").write_text("\n".join(self.genes) + "\n")
        (L.H1BANK / "perts.txt").write_text("A\nB\nC\nZ\n")
        np.save(L.H1BANK / "delta.npy", np.zeros((4, 4), np.float32))
        np.save(L.H1BANK / "b_ctrl.npy", np.zeros(4, np.float32))
        np.save(L.H1BANK / "measured_mask.npy", np.ones((4, 4), bool))
        (L.WORK / "AXIS.txt").write_text("\n".join(self.genes) + "\n")
        (L.WORK / "COHORT.tsv").write_text(
            "order\tsymbol\tr5_source_row\tbank_row\town_column_on_axis\tin_302_cohort\tin_db_joung\n"
            "0\tA\t1\t0\t1\t1\t0\n1\tB\t2\t1\t1\t0\t0\n2\tC\t3\t2\t1\t0\t0\n")
        L.SYM302.parent.mkdir(parents=True)
        L.SYM302.write_text("symbol,cohort\nA,primary\nB,secondary\n")
        (L.WORK / "R5-ROWS.tsv").write_text(
            "order\tsource_row_index\tobs_id\tsymbol\tensembl\tsign\tin_bank_cohort\n"
            "0\t1\thESC_TF_screen_A\tA\tENSGA\t1\t1\n"
            "1\t2\thESC_TF_screen_B\tB\tENSGB\t1\t1\n"
            "2\t3\thESC_TF_screen_C\tC\tENSGC\t1\t1\n")
        with h5py.File(L.DB, "r") as f:
            self.r5, _ = A._r5_identity(A._metadata(f))
        self.meta = dict(
            R5=copy.deepcopy(self.r5), file=dict(n_obs=4, n_var=4),
            cohorts=dict(nakatake_x_bank=dict(n=3, sha256=L.sha256_strings(["A", "B", "C"]))),
            axis=dict(n=4, sha256=L.sha256_strings(self.genes)))
        self.save_json(L.WORK / "meta303.json", self.meta)
        self.ka = {n: dict(pass_=True) for n in
                   ("KA-C", "KA-L1", "KA-L2", "KA-L3", "KA-L4", "KA-L5", "KA-W")}
        self.ka["summary"] = dict(all_pass=True)
        self.save_json(L.WORK / "ka303.json", self.ka)
        self.save_json(L.WORK / "test303.json",
                       dict(n=1, n_fail=0, results=[dict(pass_=True)]))
        for name in A.required_files():
            path = self.root / name
            if not path.exists():
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"synthetic dependency only\n")
        self.analysis = L.HERE / "final_analysis.py"
        self.features = L.WORK / "frozen_features.npz"
        self.analysis.write_text("# Synthetic final analysis, never executed.\n")
        np.savez(self.features, feature=np.arange(3))
        files = {name: L.sha256_file(self.root / name) for name in A.required_files()}
        files.update({A._relative(p): L.sha256_file(p) for p in (self.analysis, self.features)})
        self.rel = dict(
            schema_version=1, experiment="303", claim="C7", reserve="R5",
            dataset=L.RESERVE_DATASET, authorization="ONE_SHOT_R5_EXPRESSION",
            issued_by="coordinator", coordinator_id="synthetic-test-coordinator",
            issued_at="2026-01-01T00:00:00Z",
            manifest=dict(files=files, analysis_code=[A._relative(self.analysis)],
                          feature_files=[A._relative(self.features)], r5=copy.deepcopy(self.r5)),
            test_reports={})
        self.freeze()
        # Even an accidentally introduced production path must fail this suite.
        original_init = A._PinnedFile.__init__

        def isolated_init(pin, path):
            self.assertTrue(path.is_relative_to(self.root), f"non-fixture file open: {path}")
            original_init(pin, path)

        self.stack.enter_context(patch.object(A._PinnedFile, "__init__", isolated_init))
        self.reads = []
        self.allow_values = False
        original_getitem = h5py.Dataset.__getitem__

        def spy(dataset, key):
            if dataset.name in ("/X", "/layers/control"):
                self.assertTrue(self.allow_values, f"expression read before authorization: {dataset.name}")
                self.assertTrue(A.MARKER.is_file(), "expression read before marker")
                self.reads.append((dataset.name, key))
            return original_getitem(dataset, key)

        self.stack.enter_context(patch.object(h5py.Dataset, "__getitem__", spy))

    @staticmethod
    def save_json(path, value):
        path.write_text(json.dumps(value, allow_nan=False))

    def freeze(self):
        """Issue TEST receipts only. Does not update inputs or live row pins."""
        core = A.manifest_sha256(self.rel["manifest"])
        reports = {}
        for suite in ("access_guard303", "frozen_analysis303"):
            path = L.WORK / (suite + "-receipt.json")
            self.save_json(path, dict(schema_version=1, manifest_sha256=core,
                                      suite=suite, passed=True, tests_run=1, failures=0, errors=0))
            reports[A._relative(path)] = L.sha256_file(path)
        self.rel["test_reports"] = reports
        self.save_json(A.RELEASE, self.rel)

    def rebind(self, path):
        self.rel["manifest"]["files"][A._relative(path)] = L.sha256_file(path)
        self.freeze()

    def refused(self, pattern=None, call=None):
        with self.assertRaisesRegex(L.Refuse, pattern or "."):
            (call or A.run_access)()
        self.assertEqual(self.reads, [])

    def test_valid_preflight_is_value_free_and_marker_free(self):
        result = A.check_gate()
        self.assertEqual(result["live"]["r5"], self.r5)
        self.assertFalse(A.MARKER.exists())
        self.assertFalse(A.PRODUCT.exists())
        self.assertEqual(self.reads, [])

    def test_direct_run_without_release_refuses_before_source_open(self):
        A.RELEASE.unlink()
        with patch.object(A, "_PinnedFile", side_effect=AssertionError("opened input")):
            self.refused("no coordinator release")
            self.assertTrue(A.check_gate(False)["locked"])

    def test_missing_final_preregistration(self):
        A.PREREG_FINAL.unlink()
        with patch.object(A, "_PinnedFile", side_effect=AssertionError("opened input")):
            self.refused("final PREREGISTRATION")

    def test_empty_final_preregistration(self):
        A.PREREG_FINAL.write_bytes(b"")
        self.rebind(A.PREREG_FINAL)
        self.refused("PREREGISTRATION.md is empty")

    def test_stale_sidecar_cannot_hide_source_mutation(self):
        (L.WORK / "source.sha256").write_text(L.sha256_file(L.DB))
        with h5py.File(L.DB, "r+") as f:
            f["X"][1, 0] = 123.0
        self.refused("fresh file digest mismatch")
        self.assertFalse(A.MARKER.exists())

    def test_live_row_metadata_not_cached_meta(self):
        with h5py.File(L.DB, "r+") as f:
            f["obs/id"][1] = b"hESC_TF_screen_Q"
        self.rebind(L.DB)  # New byte hash but deliberately stale row identity.
        self.refused("live HDF5 R5 identity")

    def test_live_row_selection_not_cached_meta(self):
        with h5py.File(L.DB, "r+") as f:
            f["obs/dataset/codes"][1] = f["obs/dataset/codes"][0]
        self.rebind(L.DB)
        self.refused("live HDF5 R5 identity")

    def test_live_var_metadata_is_bound(self):
        with h5py.File(L.DB, "r+") as f:
            f["var/gene_name/categories"][0] = b"Q"
        self.rebind(L.DB)
        self.refused("live HDF5 R5 identity")

    def test_negative_categorical_code_is_not_last_category(self):
        with h5py.File(L.DB, "r+") as f:
            f["obs/dataset/codes"][1] = -1
        self.rebind(L.DB)
        self.refused("invalid categorical codes")

    def test_cohort_changed_bytes(self):
        p = L.WORK / "COHORT.tsv"
        p.write_text(p.read_text().replace("0\tA\t1", "0\tA\t2"))
        self.refused("fresh file digest mismatch")

    def test_rehashed_cohort_still_requires_live_row_mapping(self):
        p = L.WORK / "COHORT.tsv"
        p.write_text(p.read_text().replace("0\tA\t1", "0\tA\t2"))
        self.rebind(p)
        self.refused("COHORT.tsv live identity mismatch")

    def test_stale_r5_row_table_even_if_rehashed(self):
        p = L.WORK / "R5-ROWS.tsv"
        p.write_text(p.read_text().replace("0\t1", "0\t2"))
        self.rebind(p)
        self.refused("R5-ROWS.tsv live identity mismatch")

    def test_axis_order_and_live_mask_checked(self):
        p = L.WORK / "AXIS.txt"
        p.write_text("B\nA\nC\nD\n")
        self.meta["axis"]["sha256"] = L.sha256_strings(["B", "A", "C", "D"])
        self.save_json(L.WORK / "meta303.json", self.meta)
        self.rebind(p)
        self.rebind(L.WORK / "meta303.json")
        self.refused("live common measured axis")

    def test_changed_mask_even_with_updated_source_digest_refuses(self):
        with h5py.File(L.DB, "r+") as f:
            f["layers/measure_mask"][1, 0] = 0
        self.rebind(L.DB)
        self.refused("live common measured axis")

    def test_missing_mask(self):
        with h5py.File(L.DB, "r+") as f:
            del f["layers/measure_mask"]
        self.rebind(L.DB)
        self.refused("missing HDF5 object")

    def test_mask_cannot_alias_expression(self):
        with h5py.File(L.DB, "r+") as f:
            del f["layers/measure_mask"]
            f["layers/measure_mask"] = f["X"]
        self.rebind(L.DB)
        self.refused("aliases a protected expression dataset")

    def test_metadata_cannot_alias_expression(self):
        with h5py.File(L.DB, "r+") as f:
            del f["obs/id"]
            f["obs/id"] = f["X"]
        self.rebind(L.DB)
        self.refused("invalid strings")

    def test_failed_known_answer_even_if_rehashed(self):
        self.ka["KA-L1"]["pass_"] = False
        self.save_json(L.WORK / "ka303.json", self.ka)
        self.rebind(L.WORK / "ka303.json")
        self.refused("failed known answer KA-L1")

    def test_failed_saved_synthetic_tests_even_if_rehashed(self):
        p = L.WORK / "test303.json"
        self.save_json(p, dict(n=1, n_fail=0, results=[dict(pass_=1)]))
        self.rebind(p)
        self.refused("failed/incomplete frozen synthetic tests")

    def test_every_mandatory_dependency_is_required(self):
        original = copy.deepcopy(self.rel)
        for name in A.required_files():
            with self.subTest(name=name):
                self.rel = copy.deepcopy(original)
                del self.rel["manifest"]["files"][name]
                self.freeze()
                self.refused("missing mandatory files")

    def test_changed_kd_feature_analysis_and_test_files(self):
        for path in (L.H1BANK / "delta.npy", L.H1BANK / "b_ctrl.npy",
                     L.H1BANK / "perts.txt", L.H1BANK / "measured_mask.npy",
                     L.CANON, L.HELD279, L.PERSRC290, self.features, self.analysis,
                     L.HERE / "test_access_guard303.py"):
            with self.subTest(path=path):
                original = path.read_bytes()
                path.write_bytes(original + b"\nchanged")
                self.refused("fresh file digest mismatch")
                path.write_bytes(original)

    def test_release_wrong_types_claim_issuer_and_hashes(self):
        original = copy.deepcopy(self.rel)
        mutations = [
            (("schema_version",), True), (("schema_version",), 1.0),
            (("claim",), "C6"), (("reserve",), "R4"), (("issued_by",), "worker"),
            (("coordinator_id",), ""), (("coordinator_id",), 7),
            (("issued_at",), "TEST"), (("issued_at",), "2026-02-30T00:00:00Z"),
            (("issued_at",), "9999-01-01T00:00:00Z"), (("manifest",), []),
            (("test_reports",), []), (("manifest", "r5", "n_rows"), True),
            (("manifest", "r5", "row_index_sha256"), "g" * 64),
            (("manifest", "r5", "symbols_sha256"), "A" * 64),
            (("manifest", "analysis_code"), []), (("manifest", "feature_files"), []),
            (("manifest", "files", A._relative(L.DB)), None),
            (("manifest", "files", A._relative(L.DB)), 0),
            (("manifest", "files", A._relative(L.DB)), "0" * 63),
        ]
        for keys, value in mutations:
            with self.subTest(keys=keys, value=value):
                self.rel = copy.deepcopy(original)
                obj = self.rel
                for key in keys[:-1]:
                    obj = obj[key]
                obj[keys[-1]] = value
                self.save_json(A.RELEASE, self.rel)
                self.refused()

    def test_json_duplicate_extra_nonfinite_and_nonobject(self):
        valid = json.dumps(self.rel)
        for raw in ('{"claim":"C7",' + valid[1:],
                    valid[:-1] + ',"extra":1}', '{"schema_version":NaN}', "[]"):
            with self.subTest(raw=raw[:50]):
                A.RELEASE.write_text(raw)
                self.refused()

    def test_missing_stale_or_failed_receipts(self):
        original = copy.deepcopy(self.rel)
        for field, value in (("passed", False), ("passed", 1), ("tests_run", True),
                             ("tests_run", 0), ("failures", 1), ("errors", True),
                             ("manifest_sha256", "0" * 64), ("suite", "unrelated")):
            with self.subTest(field=field, value=value):
                self.rel = copy.deepcopy(original)
                self.freeze()
                name = next(iter(self.rel["test_reports"]))
                path = self.root / name
                rec = json.loads(path.read_text())
                rec[field] = value
                self.save_json(path, rec)
                self.rel["test_reports"][name] = L.sha256_file(path)
                self.save_json(A.RELEASE, self.rel)
                self.refused()
        self.rel = copy.deepcopy(original)
        self.freeze()
        self.rel["test_reports"].pop(next(iter(self.rel["test_reports"])))
        self.save_json(A.RELEASE, self.rel)
        self.refused("both required")

    def test_receipt_cycle_and_path_escape_refused(self):
        original = copy.deepcopy(self.rel)
        for name in (next(iter(self.rel["test_reports"])), A._relative(A.RELEASE),
                     "../escape", "/absolute/path", "x/../escape"):
            with self.subTest(name=name):
                self.rel = copy.deepcopy(original)
                self.rel["manifest"]["files"][name] = "0" * 64
                self.save_json(A.RELEASE, self.rel)
                self.refused()

    def test_no_reuse_of_successful_check(self):
        A.check_gate()
        A.RELEASE.unlink()
        self.refused("no coordinator release")

    def test_empty_and_dangling_marker_refuse(self):
        A.MARKER.write_bytes(b"")
        self.refused("SPENT")
        A.MARKER.unlink()
        A.MARKER.symlink_to(L.WORK / "missing")
        self.refused("SPENT")

    def test_existing_product_without_marker_refuses(self):
        A.PRODUCT.write_bytes(b"incomplete")
        self.refused("output already exists")

    def test_success_marker_fsync_before_values_and_second_call_refuses(self):
        events = []
        original_fsync = os.fsync
        original_getitem = h5py.Dataset.__getitem__

        def fsync(fd):
            mode = os.fstat(fd).st_mode
            if A.MARKER.exists():
                events.append("directory" if stat.S_ISDIR(mode) else "file")
            return original_fsync(fd)

        def read(dataset, key):
            if dataset.name in ("/X", "/layers/control"):
                self.assertIn("file", events)
                self.assertIn("directory", events)
                self.assertLess(events.index("file"), events.index("directory"))
            return original_getitem(dataset, key)

        self.allow_values = True
        with patch.object(os, "fsync", fsync), patch.object(h5py.Dataset, "__getitem__", read):
            result = A.run_access()
        self.assertTrue(result["cpm_gate_pass"])
        self.assertEqual(result["n_rows"], 3)
        with np.load(A.PRODUCT, allow_pickle=False) as out:
            np.testing.assert_array_equal(out["rows"], self.rows)
            expected = (np.log1p(5 * self.cpm[1:]) - np.log1p(12500)).astype(np.float32)
            np.testing.assert_allclose(out["delta"], expected, atol=1e-6, rtol=0)
        self.assertEqual(len(self.reads), 4)  # two row blocks, X + control per block
        self.reads.clear()
        self.refused("SPENT")

    def test_interrupted_first_slice_still_spent(self):
        original_getitem = h5py.Dataset.__getitem__

        def interrupted(dataset, key):
            if dataset.name == "/X":
                raise KeyboardInterrupt("synthetic interruption after marker")
            return original_getitem(dataset, key)

        with patch.object(h5py.Dataset, "__getitem__", interrupted):
            with self.assertRaises(KeyboardInterrupt):
                A.run_access()
        self.assertTrue(A.MARKER.exists())
        self.assertFalse(A.PRODUCT.exists())
        self.refused("SPENT")

    def test_interrupted_after_values_still_spent(self):
        self.allow_values = True
        with patch.object(L, "bank_transform", side_effect=KeyboardInterrupt("after values")):
            with self.assertRaises(KeyboardInterrupt):
                A.run_access()
        self.assertEqual(len(self.reads), 2)
        self.assertTrue(A.MARKER.exists())
        self.assertFalse(A.PRODUCT.exists())
        self.reads.clear()
        self.refused("SPENT")

    def test_marker_file_or_directory_fsync_failure_never_reads_values(self):
        original_fsync = os.fsync
        for directory in (False, True):
            with self.subTest(directory=directory):
                def failure(fd):
                    if A.MARKER.exists() and stat.S_ISDIR(os.fstat(fd).st_mode) == directory:
                        raise OSError("synthetic fsync failure")
                    return original_fsync(fd)
                with patch.object(os, "fsync", failure):
                    self.refused("fsync failure")
                self.assertTrue(A.MARKER.exists())
                self.refused("SPENT")
                A.MARKER.unlink()  # TEST fixture only, to exercise the other branch.

    def test_marker_race_is_exclusive(self):
        original = A._exclusive_json

        def competitor(path, payload):
            A.MARKER.write_bytes(b"another synthetic accessor won")
            return original(path, payload)

        with patch.object(A, "_exclusive_json", competitor):
            self.refused()
        self.assertEqual(A.MARKER.read_bytes(), b"another synthetic accessor won")

    def test_source_replaced_after_marker_before_first_value(self):
        original = A._exclusive_json

        def replacement(path, payload):
            original(path, payload)
            other = L.DB.with_suffix(".replacement")
            shutil.copyfile(L.DB, other)
            os.replace(other, L.DB)  # identical bytes, DIFFERENT inode

        with patch.object(A, "_exclusive_json", replacement):
            self.refused("file changed")
        self.assertTrue(A.MARKER.exists())
        self.refused("SPENT")

    def test_inplace_change_during_hash_detected(self):
        original = A._PinnedFile.unchanged
        counts = {}

        def mutate(pin):
            counts[pin.path] = counts.get(pin.path, 0) + 1
            if pin.path == L.DB and counts[pin.path] == 3:
                with open(L.DB, "ab") as fh:
                    fh.write(b"changed during digest")
            original(pin)

        with patch.object(A._PinnedFile, "unchanged", mutate):
            self.refused("file changed")
        self.assertFalse(A.MARKER.exists())

    def test_changed_input_after_marker_detected(self):
        original = A._exclusive_json

        def mutation(path, payload):
            original(path, payload)
            self.features.write_bytes(b"new features")

        with patch.object(A, "_exclusive_json", mutation):
            self.refused("file changed")
        self.assertTrue(A.MARKER.exists())

    def test_symlink_source_refused(self):
        other = L.DB.with_suffix(".actual")
        L.DB.rename(other)
        L.DB.symlink_to(other)
        self.refused("symlink path refused")

    def test_external_hdf5_link_refused(self):
        external = self.root / "external.h5"
        with h5py.File(external, "w") as f:
            f["X"] = np.zeros((4, 4))
        with h5py.File(L.DB, "r+") as f:
            del f["X"]
            f["X"] = h5py.ExternalLink(str(external), "X")
        self.rebind(L.DB)
        self.refused("non-local")

    def test_rehashed_loaded_code_change_refuses_restart_required(self):
        p = L.HERE / "lib303.py"
        p.write_bytes(p.read_bytes() + b"\n# changed after import\n")
        self.rebind(p)
        self.refused("loaded lib303.py differs")


class RealInputSchemaTests(unittest.TestCase):
    def test_saved_metadata_schema_and_302_overlap_without_source_or_release(self):
        """Use actual headers, not an invented fixture, with zero HDF5/source I/O."""
        with patch.object(h5py, "File", side_effect=AssertionError("no real HDF5 open")), \
                patch.object(A, "_PinnedFile", side_effect=AssertionError("no access inputs")), \
                patch.object(A, "_exclusive_json", side_effect=AssertionError("no access writes")):
            raw = REAL_SYM302.read_bytes()
            reader = csv.DictReader(io.StringIO(raw.decode("utf-8")))
            self.assertIn("symbol", reader.fieldnames)
            self.assertNotIn("tf", reader.fieldnames)
            self.assertIn("S_g", reader.fieldnames)
            rows = list(reader)
            self.assertEqual(len(rows), 175)
            self.assertTrue(all(r["S_g"].strip() for r in rows))
            primary = A._primary302(raw)  # The actual parser used by _preflight.
            r5 = list(csv.DictReader(
                io.StringIO((REAL_WORK / "R5-ROWS.tsv").read_text()), delimiter="\t"))
            cohort = list(csv.DictReader(
                io.StringIO((REAL_WORK / "COHORT.tsv").read_text()), delimiter="\t"))
            r5_symbols = {r["symbol"] for r in r5}
            self.assertEqual(len(r5_symbols & primary), 45)
            self.assertEqual(len(r5_symbols & {r["symbol"] for r in rows}), 58)
            for row in cohort:
                self.assertEqual(row["in_302_cohort"], str(int(row["symbol"] in primary)))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path,
                    help="core manifest to bind the stdout-only receipt to")
    args = ap.parse_args()
    manifest = None
    if args.manifest:
        manifest = A._json(args.manifest.read_bytes())
        # Only code files are opened here, NEVER the real source or reserve.
        for path in (REAL_HERE / "access303.py", REAL_HERE / "lib303.py",
                     Path(__file__).resolve()):
            expected = manifest["files"].get(A._relative(path))
            A._hash(expected, str(path))
            A._require(L.sha256_file(path) == expected, f"manifest does not bind tested code: {path}")
    result = unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite([
        unittest.defaultTestLoader.loadTestsFromTestCase(AccessGuardTests),
        unittest.defaultTestLoader.loadTestsFromTestCase(RealInputSchemaTests),
    ]))
    if manifest is not None:
        print(json.dumps(dict(
            schema_version=1, manifest_sha256=A.manifest_sha256(manifest),
            suite="access_guard303", passed=result.wasSuccessful(), tests_run=result.testsRun,
            failures=len(result.failures), errors=len(result.errors)), sort_keys=True))
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
