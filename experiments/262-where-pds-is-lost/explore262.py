"""262 -- first look at the per-fold diagnostics (prints only; the tables come from analyse262.py)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CANCER = FOLDS[:3]

D = {f: np.load(WORK / f"diag262-{f}.npz") for f in FOLDS}
M = {f: json.loads((WORK / f"diag262-{f}.json").read_text()) for f in FOLDS}

for f in FOLDS:
    z = D[f]
    t = z["targets"]
    p = z["cd4svd5|pds"]
    n = len(t)
    q = np.quantile(p, [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])
    dec = int(np.ceil(0.1 * n))
    order = np.argsort(p)
    bottom = order[:dec]
    print(f"\n=== {f}: n={n} mean pds {p.mean():.5f} deficit {1-p.mean():.5f}; quantiles {np.round(q, 4)}; "
          f"at 1.0: {(p == 1.0).sum()}; <0.5: {(p < 0.5).sum()}; <0.9: {(p < 0.9).sum()}")
    print("  bottom decile share of deficit:", round(float((1 - p[bottom]).sum() / (1 - p).sum()), 3),
          " bottom-quartile share:", round(float((1 - p[order[:n // 4]]).sum() / (1 - p).sum()), 3))
    print("  bottom decile:", ", ".join(f"{t[i]}({p[i]:.2f})" for i in bottom))
    # beat concentration on the truth side
    B = z["cd4svd5|beats"].astype(np.int64)
    hub = B.sum(0)            # how many (prediction, frame) pairs each truth beats
    hs = np.sort(hub)[::-1]
    print(f"  beats total {B.sum()}; top-5 truths share {hs[:5].sum()/B.sum():.3f}; top-10 {hs[:10].sum()/B.sum():.3f}; "
          f"top-30 {hs[:30].sum()/B.sum():.3f}; truths with 0 beats {(hub == 0).sum()}")
    top_hubs = np.argsort(hub)[::-1][:10]
    print("  hub truths:", ", ".join(f"{t[i]}({hub[i]})" for i in top_hubs))
    # prediction-side concentration
    rowb = B.sum(1)
    rs = np.sort(rowb)[::-1]
    print(f"  bottom-10 predictions share of beats {rs[:10].sum()/B.sum():.3f}, bottom-27 {rs[:27].sum()/B.sum():.3f}")
    # covariates vs pds
    covs = [k for k in z.files if k.startswith("cov|")]
    geo = ["cd4svd5|own", "cd4svd5|comp_max", "cd4svd5|comp_med", "cd4svd5|comp_mean", "cd4svd5|pred_nn",
           "cd4svd5|pred_mean", "cd4svd5|truth_nn", "cd4svd5|truth_mean", "cd4svd5|pred_norm", "cd4svd5|truth_norm",
           "cd4svd5|own_vs_nn", "BASE|pds", "CD4ONLY|pds", "UNIONONLY|pds", "h1union|pds", "cd4svd5|effect_pds"]
    rows = []
    for k in covs + geo:
        v = z[k].astype(np.float64)
        ok = np.isfinite(v)
        if ok.sum() < 10 or np.nanstd(v[ok]) == 0:
            rows.append((k, float("nan"), float("nan")))
            continue
        rho = spearmanr(v[ok], p[ok]).correlation
        # AUC: P(cov(bottom) > cov(rest))
        b = np.zeros(n, bool); b[bottom] = True
        vb, vr = v[b & ok], v[~b & ok]
        auc = float((vb[:, None] > vr[None, :]).mean() + 0.5 * (vb[:, None] == vr[None, :]).mean())
        rows.append((k, float(rho), auc))
    rows.sort(key=lambda r: -abs(r[1]) if np.isfinite(r[1]) else 0)
    print("  covariate ranking (spearman with pds; AUC bottom-decile-higher):")
    for k, rho, auc in rows:
        print(f"    {k:32s} rho {rho:+.3f}  auc {auc:.3f}")
    # gain/loss from h1union
    dp = z["cd4svd5|pds"] - z["h1union|pds"]
    print(f"  h1union->cd4svd5: mean {dp.mean():+.5f}; gainers {(dp > 0).sum()} losers {(dp < 0).sum()} same {(dp == 0).sum()}; "
          f"sum gains {dp[dp > 0].sum():+.3f} sum losses {dp[dp < 0].sum():+.3f}")
    o = np.argsort(dp)
    print("   top losers:", ", ".join(f"{t[i]}({dp[i]:+.3f})" for i in o[:8]))
    print("   top gainers:", ", ".join(f"{t[i]}({dp[i]:+.3f})" for i in o[::-1][:8]))
    print("  support:", M[f]["meta"]["n_whitened_hist"], "cd4:", M[f]["meta"]["cd4_count_hist"],
          "union:", M[f]["meta"]["union_count_hist"], "essential:", M[f]["meta"]["n_essential_on_panel"])

# fold agreement of the bottom decile (cancer folds)
sets = {}
for f in CANCER:
    z = D[f]; p = z["cd4svd5|pds"]; t = z["targets"]
    dec = int(np.ceil(0.1 * len(t)))
    sets[f] = set(t[np.argsort(p)[:dec]])
common = set.intersection(*[set(D[f]["targets"]) for f in CANCER])
print("\ncommon targets across cancer folds:", len(common))
for i, a in enumerate(CANCER):
    for b in CANCER[i + 1:]:
        A, Bs = sets[a] & common, sets[b] & common
        print(f"  {a} vs {b}: |A|={len(A)} |B|={len(Bs)} overlap {len(A & Bs)} jaccard {len(A & Bs)/len(A | Bs):.3f} "
              f"expected overlap {len(A)*len(Bs)/len(common):.1f}")
allthree = sets[CANCER[0]] & sets[CANCER[1]] & sets[CANCER[2]]
two = {x for x in common if sum(x in sets[f] for f in CANCER) >= 2}
print("  in all three:", sorted(allthree)); print("  in >=2:", sorted(two))
# rank correlation of per-target pds across folds
for i, a in enumerate(CANCER):
    for b in CANCER[i + 1:]:
        ta, tb = list(D[a]["targets"]), list(D[b]["targets"])
        ia = {x: j for j, x in enumerate(ta)}; ib = {x: j for j, x in enumerate(tb)}
        cm = sorted(set(ta) & set(tb))
        pa = D[a]["cd4svd5|pds"][[ia[x] for x in cm]]; pb = D[b]["cd4svd5|pds"][[ib[x] for x in cm]]
        print(f"  spearman(pds_{a}, pds_{b}) = {spearmanr(pa, pb).correlation:+.3f} on {len(cm)}")
