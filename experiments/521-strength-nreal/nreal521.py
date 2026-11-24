"""521: does 518's donor response-strength key predict a target's DE-pool size n_real in a held-out line better than
192's predictor (LOTO Pearson +0.24 on H1)? If it did, 192's per-target call-set sizing (oracle worth +0.036 overall,
unbuildable at r 0.24) would reopen. H1 truth; strength from K562-GWPS / Orion HCT116 / Orion HEK293T (518's
strength_for, noise-normalised count of strongly responding genes, never H1 data).
    uv run python experiments/521-strength-nreal/nreal521.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/518-bulk-mse-gap"))
sys.path.insert(0, str(ROOT / "experiments/513-reach-walk-anatomy"))
from h1skey518 import BANKS, strength_for  # noqa: E402
from walk513 import ALPHA, load  # noqa: E402


def rank(v):
    r = np.full(v.shape, np.nan); ok = np.isfinite(v)
    r[ok] = np.argsort(np.argsort(v[ok])) / max(ok.sum() - 1, 1)
    return r


d = load(); T = list(d["targets"]); G = d["genes"]
nreal = np.array([int(((d["rq"][k] < ALPHA) & (G != t)).sum()) for k, t in enumerate(T)])
cols = {}
for k, p in BANKS.items():
    v = strength_for(p, T)
    cols[k] = v
    ok = np.isfinite(v)
    r = np.corrcoef(v[ok], np.log1p(nreal[ok]))[0, 1]
    rs = np.corrcoef(rank(v)[ok], rank(np.log1p(nreal).astype(float))[ok])[0, 1]
    print(f"{k:14s} covers {ok.sum():2d}/85  Pearson(strength, log1p n_real) {r:+.3f}  Spearman {rs:+.3f}")
z = np.nanmean(np.vstack([rank(v) for v in cols.values()]), 0)
ok = np.isfinite(z)
print(f"{'mean rank':14s} covers {ok.sum():2d}/85  Pearson(z, log1p n_real) {np.corrcoef(z[ok], np.log1p(nreal[ok]))[0,1]:+.3f}  "
      f"Spearman {np.corrcoef(rank(z)[ok], rank(np.log1p(nreal).astype(float))[ok])[0,1]:+.3f}")
# how well does it separate the 'negligible' third (n_real < 20) from the rest?
neg = nreal < 20
auc = np.mean([(z[i] < z[j]) + 0.5 * (z[i] == z[j]) for i in np.flatnonzero(neg & ok) for j in np.flatnonzero(~neg & ok)])
print(f"AUC(strength ranks a non-negligible target above a negligible one): {auc:.3f}  (negligible n={int((neg & ok).sum())})")
for q in (0.2, 0.33, 0.5):
    lo = z <= np.nanquantile(z, q)
    print(f"  bottom {int(q*100)}% by strength: n={int((lo & ok).sum())}, of which n_real<20: {int((lo & ok & neg).sum())} "
          f"(base rate {neg[ok].mean():.2f}); median n_real {np.median(nreal[lo & ok]):.0f} vs top {np.median(nreal[~lo & ok]):.0f}")
