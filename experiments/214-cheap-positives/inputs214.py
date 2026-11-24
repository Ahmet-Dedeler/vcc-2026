"""214 — rule 11: a path and a sha256 for every operator and frozen input, never a
description.

AGENTS §11.1 rule 11.  Fifteen whitening operators now exist on disk under near-identical
names; 200 found a pair agreeing to five significant figures on the max eigenvalue and
differing only in sha, and 208 found a second pair (199's D/E/F refit) sharing the shipped
operator's name pattern, shape AND eigenvalue scale.  Production is 29,394 x 17,942 with
sha `e2f65ecf...`.

    uv run python inputs214.py
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

SHIPPED = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"

FROZEN = [
    # the two artifacts this experiment prices against and reads
    "submissions/210-cd4f034.h5ad.meta.json",
    "submissions/202-rc_callset.h5ad.meta.json",
    "submissions/a15-cd4f010.vcc.meta.json",
    "submissions/a23-massfix.h5ad.meta.json",
    # the harness and its reference
    "data/harness/h1-test.h5ad",
    "experiments/191-rc-hybrid-build/work/de191.a15.npz",
    "experiments/202-callset-rider/work/de202.rc_callset.npz",
    "experiments/169-sharpening-six-member/work/h1test-effect.npz",
    "experiments/185-rank-controlled-emission/work/a15_totals.npz",
    # the PLOLO fold caches and the challenge-panel parts (200's decomposition)
    "experiments/200-production-ablation/work/parts-k562_gwps.npz",
    "experiments/200-production-ablation/work/parts-orion_hct116.npz",
    "experiments/200-production-ablation/work/parts-orion_hek293t.npz",
    "experiments/200-production-ablation/work/parts-kolf2.npz",
    "experiments/200-production-ablation/work/parts-h1test.npz",
    "experiments/200-production-ablation/work/parts-challenge.npz",
    "experiments/200-production-ablation/work/cd4-challenge.npz",
    "experiments/200-production-ablation/work/cd4-h1test.npz",
    # 096's board solve and the mse anchors
    "experiments/096-mse-floor/work/board_solve.npz",
    "experiments/096-mse-floor/work/v11-stats.npz",
    "experiments/096-mse-floor/work/cp-proxy-ntc.npz",
    "experiments/178-composed-riders/work/eff.a15.npy",
    "experiments/196-emitted-depth/work/chal-totals-a15.npz",
    # the PerturbMatch bank itself
    "data/external/perturbmatch/delta.npy",
    "data/external/perturbmatch/genes.txt",
    "data/external/perturbmatch/perts.txt",
    # code this experiment executes but does not own
    "src/vcc2026/donor_transfer.py",
    "src/vcc2026/rc_emit.py",
    "scripts/make_submission_pm_cd4.py",
    "experiments/162-retrodicting-validator/plolo.py",
    "experiments/200-production-ablation/folds200.py",
    "experiments/200-production-ablation/h1parts200.py",
    "experiments/200-production-ablation/mse200.py",
    "experiments/196-emitted-depth/ham196.py",
    "experiments/196-emitted-depth/recon196.py",
    "experiments/185-rank-controlled-emission/harness185.py",
]


def sha256(path: Path, chunk=1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def main() -> None:
    ops = sorted({p for pat in ("data/donors/whitener*.npz",
                                "experiments/**/whitener*.npz",
                                "data/donors/**/whitener*.npz")
                  for p in ROOT.glob(pat)})
    print(f"{'path':<66}{'sha256[16]':<18}{'max eig':>10}{'n':>8}")
    rows = []
    for p in ops:
        s = sha256(p)
        try:
            with np.load(p) as z:
                ev = np.asarray(z["evals"])
                mx, n = float(ev.max()), int(ev.size)
        except Exception as e:                                   # noqa: BLE001
            mx, n = float("nan"), -1
            print(f"  (unreadable: {e})")
        mark = "  <== SHIPPED" if s == SHIPPED else ""
        rows.append(dict(path=str(p.relative_to(ROOT)), sha256=s, max_eig=mx, n=n,
                         shipped=s == SHIPPED))
        print(f"{str(p.relative_to(ROOT)):<66}{s[:16]:<18}{mx:>10.6f}{n:>8}{mark}")
    n_ship = sum(r["shipped"] for r in rows)
    print(f"\n{len(rows)} operators on disk, {len(set(r['sha256'] for r in rows))} distinct "
          f"sha256, {n_ship} matching the shipped `e2f65ecf...`")
    if n_ship != 1:
        raise SystemExit(f"ABORT (rule 11): {n_ship} files match the shipped operator sha")

    print(f"\n{'frozen input':<64}{'sha256[16]':<18}{'bytes':>16}")
    fr = []
    for rel in FROZEN:
        p = ROOT / rel
        if not p.exists():
            print(f"{rel:<64}{'MISSING':<18}")
            fr.append(dict(path=rel, missing=True))
            continue
        s = sha256(p)
        fr.append(dict(path=rel, sha256=s, bytes=p.stat().st_size))
        print(f"{rel:<64}{s[:16]:<18}{p.stat().st_size:>16,}")

    json.dump(dict(operators=rows, frozen=fr, shipped_sha256=SHIPPED),
              open(WORK / "inputs214.json", "w"), indent=1, default=float)
    print(f"\nwrote {WORK/'inputs214.json'}")


if __name__ == "__main__":
    main()
