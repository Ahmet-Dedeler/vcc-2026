#!/usr/bin/env python3
"""Read-only September 8 sweep audit; extract one payload at a time to temporary scratch.

Output is evidence, never deletion permission. Historical source files are gone, so
comparisons establish agreement with retained pre-deletion digests, not new access
to those source bytes. The 255 row has no historical content digest: report that gap.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
from pathlib import Path
import re
import tempfile

from check_submission_novelty import _hash_h5ad, sha256
from restore_prediction_h5ad import restore, storage_content_sha256

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "docs/storage-sweep-2026-09-08"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--scratch", type=Path, default=ROOT / "work")
    args = ap.parse_args()
    lossless = list(csv.DictReader((MANIFEST / "lossless.tsv").open(), delimiter="\t"))
    superseded = list(csv.DictReader((MANIFEST / "superseded.tsv").open(), delimiter="\t"))
    ledger = json.loads((ROOT / "submissions/SUBMITTED_LEDGER.json").read_text())
    report = {"started_at": dt.datetime.now(dt.timezone.utc).isoformat(),
              "manifest_sha256": {p.name: sha256(p) for p in MANIFEST.glob("*.tsv")},
              "ledger_records_before": len(ledger), "packages": [],
              "lossless_rows": len(lossless), "superseded_rows": len(superseded),
              "lossless_nominal_bytes": sum(int(r["bytes"]) for r in lossless),
              "superseded_nominal_bytes": sum(int(r["bytes"]) for r in superseded),
              "manifest_paths_still_present": [
                  r["path"] for r in lossless + superseded if (ROOT / r["path"]).exists()]}
    captures = {}
    for path in (ROOT / "data/board-entries").glob("*.json"):
        rec = json.loads(path.read_text())
        if rec.get("package"):
            captures[rec["package"]] = (path, rec)
    all_packages = sorted((ROOT / "submissions").glob("*.vcc"))
    args.scratch.mkdir(parents=True, exist_ok=True)
    for package in all_packages:
        rel = str(package.relative_to(ROOT))
        rows = [r for r in lossless if r["package"] == rel]
        digest = sha256(package)
        old = next(((h, r) for h, r in ledger.items() if r["package"] == rel), None)
        item = {"package": rel, "bytes": package.stat().st_size, "sha256": digest,
                "ledger_package_match": digest == old[0] if old else None,
                "restorable_source_paths": [r["path"] for r in rows]}
        if old or rows:
            with tempfile.TemporaryDirectory(prefix="sweep-audit-", dir=args.scratch) as scratch:
                raw = Path(scratch) / "payload.h5ad"
                restore(package, raw, raw=True)
                item["novelty_content_sha256"] = _hash_h5ad(raw)
                if old:
                    item["ledger_content_match"] = (
                        item["novelty_content_sha256"] == old[1].get("content_sha256"))
                if rows:
                    item["storage_content_sha256"] = storage_content_sha256(raw)
                    item["storage_manifest_checks"] = [
                        {"path": r["path"],
                         "match": item["storage_content_sha256"] == r["content_sha256"]
                         if re.fullmatch(r"[0-9a-f]{64}", r["content_sha256"]) else None}
                        for r in rows]
                item["payload_readable"] = True
        meta = package.with_suffix(".h5ad.meta.json")
        if meta.exists():
            item["source_h5ad_sha256"] = json.loads(meta.read_text()).get("output", {}).get("sha256")
        if rel in captures:
            path, rec = captures[rel]
            entry = rec.get("entry", {})
            item["capture"] = str(path.relative_to(ROOT))
            item["captured_entry_id"] = entry.get("entry_id", entry.get("id"))
            item["captured_status"] = entry.get("status")
            item["captured_package_sha256"] = rec.get("sha256")
        report["packages"].append(item)
        print(f"checked {rel}: package={digest[:16]} "
              f"payload={'yes' if item.get('payload_readable') else 'not requested'}", flush=True)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2) + "\n")
    report["finished_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
    report["missing_lossless_packages"] = sorted({
        r["package"] for r in lossless if not (ROOT / r["package"]).is_file()})
    report["missing_ledger_packages"] = sorted({
        r["package"] for r in ledger.values() if not (ROOT / r["package"]).is_file()})
    report["historical_digest_gaps"] = [r["path"] for r in lossless
                                       if not re.fullmatch(r"[0-9a-f]{64}", r["content_sha256"])]
    report["failed_checks"] = [
        p["package"] for p in report["packages"]
        if p.get("ledger_package_match") is False or p.get("ledger_content_match") is False
        or any(c["match"] is False for c in p.get("storage_manifest_checks", []))]
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    return int(bool(report["failed_checks"] or report["missing_lossless_packages"]
                    or report["missing_ledger_packages"]))


if __name__ == "__main__":
    raise SystemExit(main())
