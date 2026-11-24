"""Coordinator-only C7 manifest, test receipts and release issuance.

This utility NEVER calls run_access or reads expression values. The source
file is hashed; HDF5 access is obs/var metadata only. Release requires the
protocol/code to be committed and a C7 ledger entry containing the manifest
hash. Every generated file is exclusive-create to preserve the audit trail.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone
import h5py

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "experiments/303-second-overexpression-library"
sys.path.insert(0, str(HERE))
import access303 as A
import lib303 as L

CORE = L.WORK / "FROZEN-MANIFEST.json"
RECEIPTS = {
    "access_guard303": L.WORK / "FROZEN-GUARD-TESTS.json",
    "frozen_analysis303": L.WORK / "FROZEN-ANALYSIS-TESTS.json",
}


def emit(path, value):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(fd, "w") as f:
        json.dump(value, f, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())


def fresh_check(manifest):
    for rel, digest in manifest["files"].items():
        if L.sha256_file(ROOT / rel) != digest:
            raise RuntimeError(f"frozen input changed: {rel}")


def prepare():
    if A.RELEASE.exists() or A.MARKER.exists() or A.PRODUCT.exists():
        raise RuntimeError("release/access already exists; no refreeze")
    code = sorted([str(p.relative_to(ROOT)) for p in HERE.glob("*.py")] +
                  [str(Path(__file__).resolve().relative_to(ROOT))])
    features = sorted(str(p.relative_to(ROOT)) for p in
                      (L.HELD279, L.PERSRC290, L.ZNF279, L.SYM302, L.SUMS302))
    paths = sorted(set(A.required_files()) | set(code) | set(features))
    r5, _ = A._r5_identity(L.read_db_meta())
    # Exercise the actual current metadata/mask contract before creating the
    # immutable manifest. No X/control element is read by these helpers.
    with ExitStack() as stack:
        pins = {}
        for rel in A.required_files():
            p = A._PinnedFile(ROOT / rel)
            stack.callback(p.fh.close)
            pins[rel] = p
        f = stack.enter_context(h5py.File(pins[str(L.DB.relative_to(ROOT))].fh, "r"))
        meta = A._metadata(f)
        live_r5, rows = A._r5_identity(meta)
        A._preflight(f, pins, meta, live_r5, rows)
        if live_r5 != r5:
            raise RuntimeError("metadata changed during preparation")
    print("real metadata/mask preflight passed without expression", flush=True)
    manifest = dict(
        files={p: L.sha256_file(ROOT / p) for p in paths},
        analysis_code=code, feature_files=features, r5=r5)
    if manifest["files"][str(L.DB.relative_to(ROOT))] != L.DB_SHA:
        raise RuntimeError("source differs from verified 115 pin")
    emit(CORE, manifest)
    print("manifest_sha256", A.manifest_sha256(manifest), flush=True)
    print("preregistration_sha256", L.sha256_file(A.PREREG_FINAL), flush=True)


def test():
    man = json.loads(CORE.read_text())
    fresh_check(man)
    guard = subprocess.run(
        [sys.executable, str(HERE / "test_access_guard303.py"), "--manifest", str(CORE)],
        cwd=ROOT, text=True, capture_output=True)
    print(guard.stderr, end="", flush=True)
    if guard.returncode:
        print(guard.stdout)
        raise RuntimeError("guard tests failed; no receipt or release")
    receipt = json.loads(guard.stdout)
    emit(RECEIPTS["access_guard303"], receipt)
    analysis = subprocess.run(
        [sys.executable, str(HERE / "test_analysis303.py"), "--manifest", str(CORE)],
        cwd=ROOT, text=True, capture_output=True)
    print(analysis.stderr, end="", flush=True)
    if analysis.returncode:
        print(analysis.stdout)
        raise RuntimeError("analysis tests failed")
    analysis_receipt = json.loads(analysis.stdout)
    if (analysis_receipt.get("passed") is not True
            or type(analysis_receipt.get("tests_run")) is not int
            or analysis_receipt["tests_run"] <= 0):
        raise RuntimeError("analysis suite lacked a passing positive test count")
    fresh_check(man)
    emit(RECEIPTS["frozen_analysis303"], analysis_receipt)
    print("both manifest-bound test receipts saved")


def issue():
    man = json.loads(CORE.read_text())
    fresh_check(man)
    protocol = str(A.PREREG_FINAL.relative_to(ROOT))
    for path in set(man["analysis_code"]) | {protocol}:
        committed = subprocess.run(["git", "show", f"HEAD:{path}"], cwd=ROOT,
                                   capture_output=True, check=True).stdout
        if hashlib.sha256(committed).hexdigest() != man["files"][path]:
            raise RuntimeError(f"not committed at the frozen version: {path}")
    digest = A.manifest_sha256(man)
    ledger = ROOT / "docs/CONFIRMATION-LEDGER.md"
    text = ledger.read_text()
    if "## Claim C7" not in text or digest not in text:
        raise RuntimeError("C7 ledger entry/manifest digest missing")
    committed_ledger = subprocess.run(
        ["git", "show", "HEAD:docs/CONFIRMATION-LEDGER.md"],
        cwd=ROOT, capture_output=True, check=True).stdout.decode()
    if digest not in committed_ledger:
        raise RuntimeError("C7 ledger entry is not committed")
    reports = {str(p.relative_to(ROOT)): L.sha256_file(p) for p in RECEIPTS.values()}
    release = dict(
        schema_version=1, experiment="303", claim="C7", reserve="R5",
        dataset=L.RESERVE_DATASET, authorization="ONE_SHOT_R5_EXPRESSION",
        issued_by="coordinator", coordinator_id="Codex takeover coordinator",
        issued_at=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        manifest=man, test_reports=reports)
    A._schema(release)
    for path in RECEIPTS.values():
        r = json.loads(path.read_text())
        if r["manifest_sha256"] != digest or r["passed"] is not True:
            raise RuntimeError(f"invalid receipt: {path}")
    emit(A.RELEASE, release)
    # The full guard still validates all receipt fields and live row mappings.
    # A failure here is pre-access: no marker and no expression read.
    report = A.check_gate()
    print("release issued; metadata-only guard passed", report["live"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("step", choices=("prepare", "test", "issue"))
    args = ap.parse_args()
    {"prepare": prepare, "test": test, "issue": issue}[args.step]()
