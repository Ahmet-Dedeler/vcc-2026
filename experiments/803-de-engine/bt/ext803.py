#!/usr/bin/env python3
"""803 §8: 811's extra public banks' rows for the two H1 surfaces' targets (h1-test 85 + h1-trainval 126), with 811's own
reader `ext811.read_bank` (imported, never edited; ext811 has no write guard and writes nothing on import). 811's compact
banks (work/banks) carry 0-11 of these targets, because ext811 extracted only the A/B/C + L1-L6 panels.
Banks: every ext811 job EXCEPT the H1 group (h1_pan2025, h1_vcc2025: the held line, never enters) and kolf2 (not on the D/E/F
day's bank list, m811.DAY_GROUPS). The three composition donors are read by the transform from data/donors, as 811 does.
Known answers: (1) every source file's sha256 == the one ext811 recorded in 811's banks/manifest.json; (2) the same reader on
the first 40 A/B/C targets reproduces 811's compact rows for them bit for bit (float32), every bank.
Writes experiments/803-de-engine/bt/banks_h1/<bank>.npz (targets, rows float32, measured) + manifest.json.
    uv run python experiments/803-de-engine/bt/ext803.py
"""
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "experiments/811-lolo-model"))
import ext811 as X  # noqa: E402  (read-only)

OUT = HERE / "banks_h1"
SKIP = {"h1_pan2025", "h1_vcc2025", "kolf2"}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    h1 = [str(t) for t in np.load(ROOT / "experiments/323-current-family-validation/work/effects.npz")["targets"]]
    tv = [str(t) for t in np.load(ROOT / "experiments/618-h1-trainval-surface/work/tv/effects.npz")["targets"]]
    abc = [r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))][:40]
    want = set(h1) | set(tv)
    man811 = json.loads((ROOT / "experiments/811-lolo-model/work/banks/manifest.json").read_text())
    jobs = {b: [D] for b, D in X.SIMPLE.items()}
    for c in X.CD4_CONDS:
        jobs[c] = [D / c for D in X.CD4_DIRS]
    man = dict(ext811_sha256=sha(ROOT / "experiments/811-lolo-model/ext811.py"), n_wanted=len(want), banks={}, ka={})
    for b, dirs in jobs.items():
        if b in SKIP:
            continue
        rows, measured, src = X.read_bank(dirs, want | set(abc))
        # KA 1: the same source files as 811's extraction
        ref_src = man811["banks"][b]["sources"]
        same_src = all(src.get(k) == v for k, v in ref_src.items() if k in src) and set(src) == set(ref_src)
        # KA 2: A/B/C rows equal 811's compact rows
        z = np.load(ROOT / f"experiments/811-lolo-model/work/banks/{b}.npz")
        pos = {str(t): i for i, t in enumerate(z["targets"])}
        chk = [t for t in abc if t in rows and t in pos]
        same_rows = all(np.array_equal(rows[t], z["rows"][pos[t]]) for t in chk) and np.array_equal(measured, z["measured"])
        man["ka"][b] = dict(sources_equal=bool(same_src), abc_rows_checked=len(chk), abc_rows_equal=bool(same_rows))
        T = sorted(t for t in rows if t in want)
        R = np.vstack([rows[t] for t in T]) if T else np.zeros((0, len(X.GENES)), np.float32)
        dst = OUT / f"{b}.npz"
        np.savez_compressed(dst, targets=np.asarray(T), rows=R, measured=measured)
        man["banks"][b] = dict(sources=src, n_targets=len(T), h1test=len(set(T) & set(h1)), h1tv=len(set(T) & set(tv)),
                               out_sha256=sha(dst))
        print(f"{b:16s} h1test {man['banks'][b]['h1test']:3d}/85 h1tv {man['banks'][b]['h1tv']:3d}/126 | KA src {same_src} "
              f"abc {len(chk)} rows equal {same_rows} ({time.time() - t0:.0f}s)", flush=True)
    man["ka_pass"] = all(v["sources_equal"] and v["abc_rows_equal"] for v in man["ka"].values())
    (OUT / "manifest.json").write_text(json.dumps(man, indent=1))
    print("KA", "PASS" if man["ka_pass"] else "FAIL")
    return 0 if man["ka_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
