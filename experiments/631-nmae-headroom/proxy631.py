#!/usr/bin/env python3
"""631 proxy: does a head-FREE selection (called & delivered DOWN, top-N by |delivered lfc|, or |lfc| >= t) capture the
head-DOWN magnitude lever? The builder module sees only the block and the controls, not the design head. Same machinery
as decomp631 (zs base, 3 seeds, 521 pools; cross-fitted global factor m on odd/even targets)."""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import decomp631 as D  # noqa: E402
L = D.L
MS = D.MS
rows = {}
for seed in (0, 1, 2):
    S = L.load_state(seed)
    T = np.nan_to_num(S["rl"])
    z = np.load(D.W630 / f"stats-{D.TAG.format(seed=seed)}_dnk250zs.npz")
    P = np.nan_to_num(np.asarray(z["lfc"], float)); q = np.asarray(z["q"], float)
    calledD = (q < 0.05) & (P < 0)
    A = np.abs(P)
    sels = {"head-DOWN (design)": S["head"] & (P < 0)}
    for N in (100, 200, 400):
        m = np.zeros_like(calledD)
        for i in range(P.shape[0]):
            idx = np.flatnonzero(calledD[i]); idx = idx[np.argsort(-A[i, idx], kind="stable")[:N]]
            m[i, idx] = True
        sels[f"called-DOWN top{N} by |l|"] = m
    for t in (0.3, 0.5, 1.0):
        sels[f"called-DOWN |l|>={t}"] = calledD & (A >= t)
    sels["called-DOWN all"] = calledD
    for f in D.REG:
        pool = L.pools(S, f)
        base, den, ok = D.nmae_rows(P, T, pool)
        even = (np.arange(T.shape[0]) % 2) == 0
        for name, Msk in sels.items():
            dm = {mm: D.nmae_rows(np.where(Msk, mm * P, P), T, pool)[0] - base for mm in MS}
            mo = min(MS, key=lambda mm: np.nanmean(dm[mm][~even])); me = min(MS, key=lambda mm: np.nanmean(dm[mm][even]))
            cf = np.where(even, dm[mo], dm[me])
            ov = (Msk & S["head"]).sum() / max(Msk.sum(), 1)
            rows.setdefault((f, name), []).append((np.nanmean(cf), np.nanmean(dm[2.0]), mo, me, Msk.sum(1).mean(), ov))
for (f, name), v in rows.items():
    a = np.array([[x[0], x[1], x[4], x[5]] for x in v])
    print(f"f {f:<4g} {name:28s} cells/target {a[:,2].mean():6.0f} head-overlap {a[:,3].mean():.2f} | cross-fit {a[:,0].mean():+.4f} "
          f"| m=2 {a[:,1].mean():+.4f} | m chosen {sorted(set([x[2] for x in v] + [x[3] for x in v]))}")
