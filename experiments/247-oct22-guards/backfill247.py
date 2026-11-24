"""247 item D — backfill `submissions/SUBMITTED_LEDGER.json` with every package scored since
2026-09-02, from `docs/BOARD-LEDGER.md`, the `BOARD-RESULT*.md` files and `docs/LOG.md`.

Each entry: the package's sha256 (hashed from the file on disk by hashes247.json, cross-checked
against every digest the record already carries), the board entry id, the submission time as
recorded, the package path, its byte size, and a note carrying the board score.  Dates are the
RECORDED submission times, not `now`, so the ledger's `list` order is the board's order.

    uv run python experiments/247-oct22-guards/backfill247.py
"""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LEDGER = ROOT / "submissions/SUBMITTED_LEDGER.json"
HASHES = json.loads((HERE / "work/hashes247.json").read_text())

# What the record says about each scored entry since 2026-09-02 (a23, 09-02 00:32, is already in).
ENTRIES = [
    {   # docs/LOG.md 3382-3388; experiments/180-def-readiness/BOARD-RESULT.md
        "package": "submissions/abc-knee.vcc",
        "entry_id": "N3ui4hCkKx9bxExG2pxE",
        "date": "2026-09-02T23:28:00+00:00",
        "note": "abc-knee: per-context knee alpha (VCC_ALPHA=knee); overall 0.1740437 (rank 24/614, "
                "-0.000043 vs a15); falsifier fired (mse -0.00226 vs prereg +0.0050); D/E/F default "
                "stays global alpha. sha cross-checked against 242 RESULTS (f5faa1e2...)",
        "expect_prefix": "f5faa1e2",
    },
    {   # experiments/185-rank-controlled-emission/BOARD-RESULT-rc_hybrid.md
        "package": "submissions/191-rc_hybrid.vcc",
        "entry_id": "nQUJZhZb3MW2uFJHBX5e",
        "date": "2026-09-03T11:20:00+00:00",
        "note": "rc_hybrid: emitted p-order at bit-identical pseudobulk; overall 0.1834190095 "
                "(+0.0093325 vs a15, rank 19/631); reach +0.054940 scaled. Package identified by "
                "path and prep time (mtime 11:18Z, submitted 11:20Z); no other digest on record",
        "expect_prefix": None,
    },
    {   # experiments/202-callset-rider/BOARD-RESULT.md
        "package": "submissions/202-rc_callset.vcc",
        "entry_id": "31oB9yWy8CSNBDEkY4iW",
        "date": "2026-09-03T17:47:00+00:00",
        "note": "rc_callset: the call set (top 400 by |H| U top 2,200 by control CPM); overall "
                "0.1918638367 (+0.0084448 vs rc_hybrid, rank 11/643); jac positive for the first "
                "time. Date = prep mtime (published ~18:05Z)",
        "expect_prefix": None,
    },
    {   # experiments/210-cd4-build/BOARD-RESULT.md
        "package": "submissions/210-cd4f034.vcc",
        "entry_id": "QSe3JlNUYZ2BzvlcsL14",
        "date": "2026-09-04T00:00:00+00:00",
        "note": "cd4f034: VCC_CD4_TARGET_F 0.10 -> 0.34; overall 0.2049195830 (+0.0130557 vs "
                "rc_callset, rank 7); armed to fire at the slot reset",
        "expect_prefix": None,
    },
    {   # experiments/222-h1-union-build/BOARD-RESULT.md, docs/BOARD-LEDGER.md
        "package": "submissions/222-h1union-g028.vcc",
        "entry_id": "iUziqE1gewfcPKzYUDaB",
        "date": "2026-09-04T03:50:39+00:00",
        "note": "h1union: VCC_H1_UNION_G=0.28 (H1 pair); overall 0.2118310 (+0.0069114 vs cd4f034, "
                "rank 3/656); THE CURRENT INCUMBENT. sha f13ce40b... recorded in BOARD-LEDGER",
        "expect_prefix": "f13ce40b",
    },
]


def main() -> None:
    ledger = json.loads(LEDGER.read_text())
    before = len(ledger)
    added = []
    for e in ENTRIES:
        h = HASHES[e["package"]]
        if e["expect_prefix"] and not h["sha256"].startswith(e["expect_prefix"]):
            raise SystemExit(f"{e['package']}: hashed {h['sha256'][:16]} but the record says "
                             f"{e['expect_prefix']}...; refusing to backfill a mismatch")
        if h["sha256"] in ledger:
            print(f"already present: {e['package']} {h['sha256'][:16]}")
            continue
        ledger[h["sha256"]] = {
            "package": e["package"], "entry_id": e["entry_id"], "date": e["date"],
            "note": e["note"], "bytes": h["bytes"],
            "backfilled_by": "experiments/247-oct22-guards/backfill247.py (2026-09-04)",
        }
        added.append((e["package"], h["sha256"]))
    LEDGER.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
    print(f"ledger {before} -> {len(ledger)} entries; added {len(added)}:")
    for pkg, sha in added:
        print(f"  {sha[:16]}  {pkg}")


if __name__ == "__main__":
    main()
