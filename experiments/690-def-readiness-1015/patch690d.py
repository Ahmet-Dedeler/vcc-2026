#!/usr/bin/env python3
"""690d — capture_board_entry.py must accept a FINAL-phase entry whose members came back without a public rank.

UNAPPLIED until the coordinator applies it (before Oct 22). Why: `validate()` requires a numeric `rank` on every published
entry and the three stamp fields (panel_id, anchor_version, partition). On the dark final board a final entry's members
may come back through `vcc status` (the CLI withholds nothing; 680's NEWS §1.2 found a gated final-week team view), while
its public rank stays null. Such a record would fail closed, and `arm_submit_after.sh` then REFUSES the day's second slot
("could not capture prior … NOT submitting"): every second final-phase slot would be lost. 531 F8 covered only the
all-null "members hidden" case.

The change, for FINAL-phase entries only (`is_final` true, or a partition other than val; the same test members_hidden
uses): a missing rank is not a problem, and a missing stamp field is recorded as a warning in the wrapper instead of a
failure. Every member must still be numeric (a partial record still fails closed). Validation (A/B/C) entries keep every
check unchanged.

    uv run python experiments/690-def-readiness-1015/patch690d.py --root <repo> [--check]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

MARK = "690d: final-phase entries"
FILE = Path("scripts/capture_board_entry.py")
EDITS = [
    ("helper",
     '''def members_hidden(entry: dict) -> bool:
''',
     '''def is_final_phase(entry: dict) -> bool:
    """690d: final-phase entries: is_final true, or a partition other than val (members_hidden's own test)."""
    return entry.get("is_final") is True or (entry.get("partition") not in (None, "", "val"))


def final_warnings(entry: dict) -> list[str]:
    """690d: what a FINAL-phase record lacks without failing: its public rank and any stamp field."""
    if not is_final_phase(entry):
        return []
    w = [f"{k} missing" for k in STAMP if not entry.get(k)]
    if entry.get("status") == "published" and not _num(entry.get("rank")):
        w.append("rank missing (dark final board)")
    return w


def members_hidden(entry: dict) -> bool:
'''),
    ("stamp",
     '''        for k in STAMP:
            if not entry.get(k):
                problems.append(f"{k} missing (scores are not comparable across panels)")
        if status == "published" and not _num(entry.get("rank")):
            problems.append("rank missing")
''',
     '''        # 690d: final-phase entries: a missing stamp or public rank is a warning (final_warnings), not a failure
        for k in STAMP:
            if not entry.get(k) and not is_final_phase(entry):
                problems.append(f"{k} missing (scores are not comparable across panels)")
        if status == "published" and not _num(entry.get("rank")) and not is_final_phase(entry):
            problems.append("rank missing")
'''),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    p = a.root.resolve() / FILE
    t = p.read_text()
    if MARK in t:
        print("REFUSED: 690d is already applied", file=sys.stderr)
        return 1
    for name, old, new in EDITS:
        if t.count(old) != 1:
            print(f"REFUSED: anchor {name} found {t.count(old)} times (want 1)", file=sys.stderr)
            return 1
        t = t.replace(old, new)
    # record the warnings in the wrapper the script writes (best effort: only if the wrapper dict literal is found)
    anchor = '"members_hidden": members_hidden(entry),'
    if t.count(anchor) == 1:
        t = t.replace(anchor, anchor + ' "final_warnings": final_warnings(entry),')
    if a.check:
        print("CHECK OK: 2 anchors (+ the wrapper's final_warnings field); nothing written")
        return 0
    p.write_text(t)
    print(f"APPLIED 690d to {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
