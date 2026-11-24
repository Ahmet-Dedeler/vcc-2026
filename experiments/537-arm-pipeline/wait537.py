#!/usr/bin/env python3
"""537 — wait for board captures: poll data/board-entries/*.json (read-only) until an entry whose `note` starts with each
requested arm label exists, or until the deadline. Prints `<label> <entry id> <score_avg> <note>` for every match.

usage: uv run python experiments/537-arm-pipeline/wait537.py <deadline YYYY-MM-DDTHH:MMZ> <label> [<label> ...]
  e.g. wait537.py 2026-10-01T06:00Z 460dr 446bj2zcp
Exit 0 when all labels matched, 2 at the deadline (partial matches printed). Polls every 5 min.
"""
from __future__ import annotations

import datetime as dt
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BE = ROOT / "data/board-entries"


def scan(labels):
    found = {}
    for f in BE.glob("*.json"):
        try:
            d = json.loads(f.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        note = str(d.get("note") or "")
        e = d.get("entry", d)
        for lab in labels:
            if note.split(" ")[0] == lab:
                found[lab] = (e.get("entry_id", f.stem), e.get("score_avg"), note)
    return found


def main():
    deadline = dt.datetime.strptime(sys.argv[1], "%Y-%m-%dT%H:%MZ").replace(tzinfo=dt.timezone.utc)
    labels = sys.argv[2:]
    while True:
        found = scan(labels)
        if len(found) == len(labels) or dt.datetime.now(dt.timezone.utc) >= deadline:
            for lab in labels:
                print(lab, *found.get(lab, ("MISSING",)), flush=True)
            sys.exit(0 if len(found) == len(labels) else 2)
        time.sleep(300)


if __name__ == "__main__":
    main()
