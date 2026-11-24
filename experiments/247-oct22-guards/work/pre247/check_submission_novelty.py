#!/usr/bin/env python3
"""Refuse to submit a package that has already been scored.

STANDING RULE (Ahmet, 2026-09-01; AGENTS.md section 6, docs/LOG.md):
    Never spend a submission slot on leaderboard rank, and never re-submit an already-scored
    package.  We get 2 slots/day and they are the only measurement of the real objective.
    Leaderboard rank is informational -- only the last entry after the Oct 22 test release is
    prize-eligible.  The goal is FIRST PLACE, not a better end-of-day row.

This exists because the rule was known and articulated in this repo hours before it was broken:
on 2026-08-31 a19 regressed, the displayed rank fell 12 -> 26, and the day's second slot went on
re-uploading a15 byte-identical -- buying rank and zero information.  Prose did not stop it, so
this is the version that can.

    # before `vcc submit`
    uv run python scripts/check_submission_novelty.py check submissions/a21-resample-x014.vcc

    # after a successful submit
    uv run python scripts/check_submission_novelty.py record submissions/a21-resample-x014.vcc \
        --entry-id QGpNVMwzzdSlhUoQtRdp --note "resampling emitter x=0.14"

`check` exits 1 if the package's sha256 is already in the ledger.  It is deliberately impossible
to pass with a flag -- if you genuinely need to re-submit identical bytes, delete the ledger
entry by hand and write down why in docs/LOG.md.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import sys
from pathlib import Path

LEDGER = Path(__file__).resolve().parent.parent / "submissions" / "SUBMITTED_LEDGER.json"


def sha256(path: Path, bufsize: int = 1 << 24) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(bufsize):
            h.update(chunk)
    return h.hexdigest()


def load() -> dict:
    if not LEDGER.exists():
        return {}
    return json.loads(LEDGER.read_text())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check", help="exit 1 if this package was already submitted")
    c.add_argument("package")
    r = sub.add_parser("record", help="record a package as submitted")
    r.add_argument("package")
    r.add_argument("--entry-id", required=True)
    r.add_argument("--note", default="")
    sub.add_parser("list", help="show the ledger")
    a = ap.parse_args()

    ledger = load()
    if a.cmd == "list":
        for h, rec in sorted(ledger.items(), key=lambda kv: kv[1].get("date", "")):
            print(f"{h[:16]}  {rec.get('date','?'):26} {rec.get('entry_id','?'):22} "
                  f"{rec.get('package','?')}  {rec.get('note','')}")
        return 0

    pkg = Path(a.package)
    if not pkg.exists():
        print(f"ERROR: no such package: {pkg}", file=sys.stderr)
        return 2
    print(f"hashing {pkg} ({pkg.stat().st_size / 1e9:.2f} GB)...", file=sys.stderr)
    digest = sha256(pkg)

    if a.cmd == "check":
        prev = ledger.get(digest)
        if prev is None:
            print(f"OK: novel package, sha256 {digest[:16]}...")
            return 0
        print("=" * 78, file=sys.stderr)
        print("REFUSED: this package has ALREADY BEEN SCORED. Do not spend a slot on it.",
              file=sys.stderr)
        print(f"  sha256   {digest}", file=sys.stderr)
        print(f"  entry    {prev.get('entry_id')}   submitted {prev.get('date')}", file=sys.stderr)
        print(f"  package  {prev.get('package')}   {prev.get('note','')}", file=sys.stderr)
        print("", file=sys.stderr)
        print("  A slot must buy information. An already-scored package provably cannot --",
              file=sys.stderr)
        print("  we know its result before we spend the slot. If a bad arm is showing on the",
              file=sys.stderr)
        print("  board, LEAVE IT: rank is informational, only the last entry after Oct 22 is",
              file=sys.stderr)
        print("  prize-eligible. See AGENTS.md section 6.", file=sys.stderr)
        print("=" * 78, file=sys.stderr)
        return 1

    ledger[digest] = {
        "package": str(pkg),
        "entry_id": a.entry_id,
        "date": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "note": a.note,
        "bytes": pkg.stat().st_size,
    }
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    LEDGER.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
    print(f"recorded {digest[:16]}... -> {a.entry_id}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
