"""295 item 4b — the PAIRED file-block bootstrap of C_held against each baseline, and of
gen_pc20 against breadth.  Reads `reserve295.csv` only (no reserve file is opened again).

`block_boot` is copied from `experiments/291-f13-independent-surfaces/lib291.py:166-195`.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

WORK = Path("~/Code/virtual-cell-challenge-2026/experiments/295-reserve-r3-claim/work")
SEED = 20260906
N_BOOT = 2000


def spearman(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return float("nan")
    return float(stats.spearmanr(a[m], b[m]).statistic)


def block_boot(y, preds: dict, groups, rng, n: int = N_BOOT) -> dict:
    groups = np.asarray(groups)
    uniq = np.array(sorted(set(groups.tolist())))
    idx_of = {g: np.flatnonzero(groups == g) for g in uniq}
    keys = list(preds)
    acc = {k: [] for k in keys}
    pair = {f"{a}-{b}": [] for a in keys for b in keys if a != b}
    for _ in range(n):
        pick = rng.integers(0, len(uniq), size=len(uniq))
        sel = np.concatenate([idx_of[uniq[p]] for p in pick])
        yy = y[sel]
        cur = {k: spearman(yy, preds[k][sel]) for k in keys}
        for k in keys:
            acc[k].append(cur[k])
        for a in keys:
            for b in keys:
                if a != b:
                    pair[f"{a}-{b}"].append(cur[a] - cur[b])
    out = {k: dict(point=spearman(y, preds[k]),
                   lo=float(np.nanpercentile(acc[k], 2.5)),
                   hi=float(np.nanpercentile(acc[k], 97.5)))
           for k in keys}
    out["_paired"] = {k: dict(mean=float(np.nanmean(v)),
                              lo=float(np.nanpercentile(v, 2.5)),
                              hi=float(np.nanpercentile(v, 97.5)),
                              p_le_0=float(np.mean(np.asarray(v) <= 0)))
                      for k, v in pair.items()}
    return out


def main() -> None:
    df = pd.read_csv(WORK / "reserve295.csv")
    P = df[df["primary"]].reset_index(drop=True)
    y = P["A"].to_numpy(float)
    preds = {k: P[k].to_numpy(float) for k in
             ["C_held", "M_ann", "mag", "r_split", "de", "r_bar", "C_held_excl", "gen_pc20",
              "breadth"]}
    rng = np.random.default_rng(SEED)
    out = block_boot(y, preds, P["file"].to_numpy(), rng)
    keep = {k: v for k, v in out["_paired"].items()
            if k.startswith("C_held-") or k in ("gen_pc20-breadth", "breadth-gen_pc20",
                                                "gen_pc20-M_ann")}
    for k, v in sorted(keep.items()):
        print(f"  {k:<24} delta {v['mean']:+.4f} [{v['lo']:+.4f}, {v['hi']:+.4f}]  "
              f"P(delta<=0) {v['p_le_0']:.3f}")
    (WORK / "paired295.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
