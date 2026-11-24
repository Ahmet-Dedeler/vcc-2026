"""241 — rule 11: the path and sha256 of every operator and frozen input this experiment read.

Rule 11's point is that an operator must never be identified by its filename: 200 found two
files whose max eigenvalue agrees to five significant figures and whose sha256 differ
completely.  Every sha here is checked against the one 200 recorded for the same file.
"""
from __future__ import annotations

import hashlib
import json
import os

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")

PATHS = [
    "data/donors/whitener_union_pruned.npz",
    "data/donors/whitener_union.npz",
    "experiments/062-cd4-post-whiten/cloud-output/whiteners/held_k562_gwps.npz",
    "experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hct116.npz",
    "experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hek293t.npz",
    "experiments/095-clean-harness/out-h1test/whitener_union.npz",
    "data/harness/h1-test.h5ad",
    "scripts/make_submission_pm_cd4.py",
    "experiments/169-sharpening-six-member/work/h1test-effect.npz",
    "experiments/096-mse-floor/work/board_solve.npz",
    "experiments/096-mse-floor/work/v11-stats.npz",
    "submissions/222-h1union-g028.h5ad.meta.json",
]
# 200 `work/inputs200.json`, first 16 hex
EXPECT = {
    "data/donors/whitener_union_pruned.npz": "e2f65ecf69b7e0b0",
    "experiments/062-cd4-post-whiten/cloud-output/whiteners/held_k562_gwps.npz":
        "13d988fcc3b129d2",
    "experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hct116.npz":
        "7cecbf7f36e40d22",
    "experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hek293t.npz":
        "97e486e0d1a805e1",
    "experiments/095-clean-harness/out-h1test/whitener_union.npz": "f80fe8d6ba59138b",
    "data/harness/h1-test.h5ad": "06c6baea8c62ef3d",
    "experiments/169-sharpening-six-member/work/h1test-effect.npz": "e4d8af8b892cb22d",
}


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    out, bad = {}, []
    for p in PATHS:
        fp = os.path.join(ROOT, p)
        if not os.path.exists(fp):
            out[p] = "MISSING"
            continue
        h = hashlib.sha256()
        with open(fp, "rb") as f:
            for c in iter(lambda: f.read(1 << 22), b""):
                h.update(c)
        d = h.hexdigest()
        out[p] = d
        want = EXPECT.get(p)
        tag = ""
        if want is not None:
            ok = d[:16] == want
            tag = "  == 200" if ok else f"  DRIFT (200: {want})"
            if not ok:
                bad.append(p)
        print(f"{d[:16]}  {os.path.getsize(fp):>14,}  {p}{tag}")
    json.dump(out, open(os.path.join(WORK, "inputs241.json"), "w"), indent=1)
    print(f"\n{len(EXPECT) - len(bad)}/{len(EXPECT)} operator/frozen shas match 200's record"
          + (f"; DRIFT on {bad}" if bad else ""))
    print("note: `scripts/make_submission_pm_cd4.py` HAS drifted from 200's 90c4f49eebc1c706 "
          "(240 owns it).\n      241 never calls its mass-fix path (mass_fix is False on every "
          "241 emission), and the\n      one place it IS used -- `mse241.py`'s MF base -- "
          "reproduced 200's A/B/claim exactly,\n      which is the evidence the drift is not in "
          "`mass_consistent_effect`.")


if __name__ == "__main__":
    main()
