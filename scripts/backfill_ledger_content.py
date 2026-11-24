#!/usr/bin/env python3
"""Backfill `content_sha256` (the prediction's content digest, computed INSIDE the container) on every
`submissions/SUBMITTED_LEDGER.json` record whose package is on disk and which lacks it (256c, 2026-09-05).

    uv run python scripts/backfill_ledger_content.py [--ledger PATH] [--only DIGESTPREFIX]

One package at a time (~1-2 min and ~13 GB of scratch each); the ledger is rewritten atomically after each
record so an interrupted run loses nothing. Records whose package is missing are left as they are and listed.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_submission_novelty as N  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ledger", type=Path, default=N.LEDGER)
    ap.add_argument("--only", default="", help="only records whose .vcc digest starts with this")
    a = ap.parse_args()
    ledger = json.loads(a.ledger.read_text())
    todo = [(h, r) for h, r in ledger.items() if not r.get("content_sha256") and (not a.only or h.startswith(a.only))]
    print(f"{len(todo)} records without content_sha256", flush=True)
    for h, rec in todo:
        pkg = Path(rec["package"])
        if not pkg.is_absolute():
            pkg = N.LEDGER.parent.parent / pkg
        if not pkg.exists():
            print(f"{h[:16]}  {rec.get('entry_id')}  {pkg.name}: package not on disk, skipped", flush=True)
            continue
        t0 = _dt.datetime.now()
        dig, src = N.content_sha256(pkg)
        dt = (_dt.datetime.now() - t0).total_seconds()
        if dig is None:
            print(f"{h[:16]}  {rec.get('entry_id')}  {pkg.name}: FAILED ({src})", flush=True)
            continue
        ledger = json.loads(a.ledger.read_text())  # re-read: never clobber a concurrent `record`
        ledger[h]["content_sha256"] = dig
        ledger[h]["content_sha256_source"] = f"backfilled {_dt.datetime.now(_dt.timezone.utc).isoformat(timespec='seconds')}: {src}"
        tmp = a.ledger.with_suffix(f".{os.getpid()}.tmp")
        tmp.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
        os.replace(tmp, a.ledger)
        print(f"{h[:16]}  {rec.get('entry_id')}  {pkg.name}: content {dig[:16]}... ({dt:.0f} s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
