#!/usr/bin/env python3
"""Print the next unused leaderboard model name from docs/MODEL-NAMES.md.

A name counts as used if it appears (case-insensitive) in:
  - submissions/SUBMITTED_LEDGER.json (any `note`),
  - data/board-entries/*.json (captured board records),
  - experiments/*/chain*.sh or arm*.sh (names baked into armed chains),
  - docs/model-names-claimed.tsv (explicit claims).

`--claim ARM` appends the chosen name to the claims file so a second chain armed the same day
cannot take it. The reserved final-entry name is never handed out.
"""

import argparse
import datetime as dt
import glob
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NAMES_MD = ROOT / "docs" / "MODEL-NAMES.md"
CLAIMS = ROOT / "docs" / "model-names-claimed.tsv"


def queue() -> list[str]:
    """Numbered entries under '## Queue', in order."""
    text = NAMES_MD.read_text()
    body = text.split("## Queue", 1)[1]
    return [m.group(1).strip() for m in re.finditer(r"^\d+\.\s+(.+)$", body, re.M)]


def used_haystack() -> str:
    """One lowercase blob of every place a name could already be used."""
    parts = []
    ledger = ROOT / "submissions" / "SUBMITTED_LEDGER.json"
    if ledger.exists():
        parts += [str(v.get("note", "")) for v in json.loads(ledger.read_text()).values()]
    for f in glob.glob(str(ROOT / "data" / "board-entries" / "*.json")):
        parts.append(Path(f).read_text())
    for pat in ("experiments/*/chain*.sh", "experiments/*/arm*.sh"):
        for f in glob.glob(str(ROOT / pat)):
            parts.append(Path(f).read_text(errors="ignore"))
    if CLAIMS.exists():
        parts.append(CLAIMS.read_text())
    return "\n".join(parts).lower()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--claim", metavar="ARM", help="reserve the name for this arm id")
    ap.add_argument("--list", action="store_true", help="show every queue name with its status")
    a = ap.parse_args()

    hay = used_haystack()
    names = queue()
    if a.list:
        for i, n in enumerate(names, 1):
            print(f"{i:3d}  {'used' if n.lower() in hay else 'free'}  {n}")
        return 0
    free = [n for n in names if n.lower() not in hay]
    if not free:
        print("queue exhausted: append screened names to docs/MODEL-NAMES.md", file=sys.stderr)
        return 1
    name = free[0]
    if a.claim:
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%MZ")
        new = not CLAIMS.exists()
        with CLAIMS.open("a") as fh:
            if new:
                fh.write("claimed_utc\tarm\tname\n")
            fh.write(f"{stamp}\t{a.claim}\t{name}\n")
    print(name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
