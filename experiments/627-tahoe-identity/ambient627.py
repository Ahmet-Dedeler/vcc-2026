#!/usr/bin/env python3
"""627b — is Tahoe's identity erosion an AMBIENT-RNA artifact of its pooled (MOSAIC, 50 lines per well) design?

Ambient soup = the mix of all lines' RNA. More soup in treated cells (death, lysis, less RNA per cell) pulls every line
toward the pooled mean, which mimics "identity genes go down". But soup has a fingerprint biology lacks: genes ABSENT in
line L (control CPM < 5, other plates) yet expressed elsewhere (mean CPM of the other lines > 50, i.e. >= 10x) gain reads in L. True
identity erosion predicts ~0 there, because cells do not start expressing other lineages' genes.
Per (d, L), with leave-plate-out basal as in tahoe627:
  r = delta - cross-line mean (line-specific deviation);
  absent_up   = mean r over L's absent-but-elsewhere-expressed genes (ambient predicts > 0, growing with strength)
  ident_down  = mean r over L's top-decile identity genes (both models predict < 0)
Read by response-strength quintile. If absent_up ~ 0 while ident_down < 0, the erosion is not ambient.
Analysis only; < 3 GB.
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W184 = ROOT / "experiments/184-tahoe-transfer/work"


def main() -> None:
    z = np.load(W184 / "tensor.npz", allow_pickle=True)
    delta = z["delta"].astype(np.float64)
    basal = z["basal_plate"].astype(np.float64)
    valid = z["ctrl_cells"] > 0
    plate_of = np.load(W184 / "plate_of.npy", allow_pickle=True).astype(np.int64)
    T, L, G = delta.shape
    P = basal.shape[1]
    cache = {}

    def prof(mask):
        if mask not in cache:
            num = np.zeros((L, G)); den = np.zeros((L, 1))
            for p in range(P):
                if (mask >> p) & 1:
                    continue
                w = valid[:, p][:, None].astype(float)
                num += np.nan_to_num(basal[:, p, :]) * w; den += w
            cpm = 1e6 * num / np.maximum(den, 1)
            lg = np.log1p(cpm)
            other_mean_cpm = np.array([np.delete(cpm, l, axis=0).mean(0) for l in range(L)])
            I = lg - np.array([np.delete(lg, l, axis=0).mean(0) for l in range(L)])
            cache[mask] = (cpm, other_mean_cpm, I)
        return cache[mask]

    rows = []
    for d in range(T):
        dd = delta[d]
        ok = np.isfinite(dd).all(1) & (np.abs(dd).sum(1) > 0)
        if ok.sum() < 20:
            continue
        r = dd - np.nanmean(dd[ok], 0)
        for l in np.flatnonzero(ok):
            cpm, oth, I = prof(int(plate_of[d, l]))
            absent = (cpm[l] < 5) & (oth[l] > 50)          # >= 10x elsewhere; < 1 CPM leaves ~0 genes on the 2,000-gene axis
            expr = cpm[l] > 5
            if absent.sum() < 5 or expr.sum() < 100:
                continue
            Ie = np.where(expr, I[l], np.nan)
            top = expr & (Ie >= np.nanquantile(Ie, 0.9))
            bot = expr & (Ie <= np.nanquantile(Ie, 0.1))
            rows.append((float(np.linalg.norm(dd[l])), float(r[l][absent].mean()), float(r[l][top].mean()),
                         float(r[l][bot].mean()), int(absent.sum())))
    a = np.array(rows)
    q = np.quantile(a[:, 0], [0.2, 0.4, 0.6, 0.8])
    out = {"n_pairs": len(a), "median_absent_genes": float(np.median(a[:, 4])), "by_strength": []}
    print(f"pairs {len(a)}; median absent-but-elsewhere-expressed genes per line {np.median(a[:, 4]):.0f}")
    print("strength Q | absent_up (ambient fingerprint) | ident_down (top-decile I) | low-I decile (expressed)")
    for k in range(5):
        s = np.digitize(a[:, 0], q) == k
        row = dict(q=k + 1, n=int(s.sum()), absent_up=float(a[s, 1].mean()), absent_up_se=float(a[s, 1].std() / np.sqrt(s.sum())),
                   ident_down=float(a[s, 2].mean()), low_I=float(a[s, 3].mean()))
        out["by_strength"].append(row)
        print(f"  Q{k + 1}  n {row['n']:6d} | {row['absent_up']:+.5f} ± {row['absent_up_se']:.5f} | {row['ident_down']:+.5f} | {row['low_I']:+.5f}")
    out["overall"] = dict(absent_up=float(a[:, 1].mean()), ident_down=float(a[:, 2].mean()), low_I=float(a[:, 3].mean()))
    print("overall:", {k: round(v, 5) for k, v in out["overall"].items()})
    (HERE / "work/ambient627.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
