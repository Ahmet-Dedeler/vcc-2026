#!/usr/bin/env python3
"""690c — keep the pinned builder a def_day run built with. UNAPPLIED until the coordinator applies it (after 690a).

Why: a D/E/F build's meta records `hashes.implementation["scripts/make_submission_pm_cd4.py"]` = the builder AFTER the
day's pins (pin-artifacts, then pin-lodo; e.g. 524dr2: d5cec53d). def_day.sh's exit trap restores the HEAD builder, and the
pin records keep only the PRE-pin copy, so the recorded builder exists nowhere. The runbook's old answer (run with
DEF_DAY_KEEP_PINS=1 and commit the pinned builder to scripts/) leaves HEAD pinned, and every later def_day run then
refuses or re-pins on top. 690c instead copies the pinned builder to <work>/builder.pinned.py inside the trap, before the
restore, and logs its sha. Committing that file with the round's records makes the meta reproducible from git while HEAD
stays clean.

One file, scripts/def_day.sh, one anchored insertion at the top of restore(). It runs only when a pin record exists (a
run that never pinned writes nothing) and never changes the trap's exit code.

    uv run python experiments/690-def-readiness-1015/patch690c.py --root <repo> [--check]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

MARK = "690c: pinned builder kept"
DAY = Path("scripts/def_day.sh")
OLD = '''restore(){
  local rc=$?
'''
NEW = '''restore(){
  local rc=$?
  # 690c: pinned builder kept -- the builder the day's builds ran on (after every pin; the final meta's
  # hashes.implementation sha), saved before the restore so the round can commit it next to its records
  if [ -f "$W/pin-artifacts.json" ] && cp scripts/make_submission_pm_cd4.py "$W/builder.pinned.py" 2>/dev/null; then
    log "690c: pinned builder saved as $W/builder.pinned.py ($(shasum -a 256 "$W/builder.pinned.py" | cut -c1-16))"
  fi
'''


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    p = a.root.resolve() / DAY
    t = p.read_text()
    if MARK in t:
        print("REFUSED: 690c is already applied", file=sys.stderr)
        return 1
    if t.count(OLD) != 1:
        print(f"REFUSED: restore() anchor found {t.count(OLD)} times (want 1)", file=sys.stderr)
        return 1
    if a.check:
        print("CHECK OK: 1 def_day.sh anchor; nothing written")
        return 0
    p.write_text(t.replace(OLD, NEW))
    print(f"APPLIED 690c to {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
