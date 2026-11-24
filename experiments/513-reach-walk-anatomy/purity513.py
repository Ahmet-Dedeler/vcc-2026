"""513e: is there an identifiable (target, gene) subset whose sign is >= 90 % right on H1's pool (reference-significant genes)?
Sources (all truth-free, H1-excluded):
  donors  same-target rows of k562_gwps, orion_hct116, orion_hek293t, rpe1 (log2 pseudobulk ratios, common480 loader);
          agreement = number of covering donors with the same sign; consensus = mean log2 ratio
  clean4  the CLEAN4 programme P_down (target-independent)
  ours    our emitted lfc sign (435 production)
Accuracy is measured on pool genes (reference q < 0.05, finite non-zero lfc, target excluded), per pool band.
    uv run python experiments/513-reach-walk-anatomy/purity513.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "480-nmae-instrument"))
from common480 import load_donor_rows  # noqa: E402
from walk513 import ALPHA, load  # noqa: E402

d = load(); G = d["genes"]; T = d["targets"]
ROOT = Path(__file__).resolve().parents[2]
z = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
pidx = {g: i for i, g in enumerate(z["genes"].astype(str))}
Pd = np.array([z["P_down"][pidx[g]] if g in pidx else np.nan for g in G])
donors = ["k562_gwps", "orion_hct116", "orion_hek293t", "rpe1"]
rows = {}
for dn in donors:
    L2, has, okc, _ = load_donor_rows(dn, T, G, pc_cpm=1.0)
    rows[dn] = L2
    print(f"{dn}: rows for {has.sum()}/{len(T)} targets, genes on axis {okc.sum()}")
S = np.stack([np.sign(rows[dn]) for dn in donors])           # donors x targets x genes (NaN where absent)
V = np.stack([rows[dn] for dn in donors])
ncov = np.isfinite(S).sum(0)
nneg = (S < 0).sum(0); npos = (S > 0).sum(0)
cons = np.nanmean(V, 0)
pool = np.zeros_like(d["rl"], bool)
for k, t in enumerate(T):
    pool[k] = (d["rq"][k] < ALPHA) & np.isfinite(d["rl"][k]) & (d["rl"][k] != 0) & (G != t)
nr = pool.sum(1)
band = np.where(nr < 20, 0, np.where(nr < 100, 1, 2))
truth = np.sign(d["rl"])
def acc(mask, pred_sign):
    m = mask & pool & np.isfinite(pred_sign) & (pred_sign != 0)
    n = int(m.sum())
    return (float((pred_sign[m] == truth[m]).mean()) if n else float("nan")), n
print("\n== donor-consensus sign on pool genes, by agreement and |consensus| (all bands; then per band)")
for need in (1, 2, 3, 4):
    for thr in (0.0, 0.5, 1.0, 2.0):
        agree = ((nneg >= need) & (npos == 0)) | ((npos >= need) & (nneg == 0))
        m = agree & (np.abs(cons) >= thr)
        a, n = acc(m, np.sign(cons))
        per = [acc(m & (band[:, None] == b), np.sign(cons)) for b in range(3)]
        print(f"  >= {need} donors, all agree, |cons| >= {thr}: acc {a:.3f} n {n:6d}   " +
              "  ".join(f"band{b} {p[0]:.3f} (n {p[1]})" for b, p in enumerate(per)))
print("\n== CLEAN4 programme sign on pool genes by confidence |P_down - 0.5|")
csign = np.where(Pd > 0.5, -1.0, np.where(Pd < 0.5, 1.0, np.nan))[None, :] * np.ones((len(T), 1))
for lo, hi in [(0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5)]:
    m = (np.abs(Pd - 0.5) >= lo)[None, :] & (np.abs(Pd - 0.5) < hi)[None, :]
    a, n = acc(m, csign)
    per = [acc(m & (band[:, None] == b), csign) for b in range(3)]
    print(f"  conf [{lo},{hi}): acc {a:.3f} n {n:6d}   " + "  ".join(f"band{b} {p[0]:.3f} (n {p[1]})" for b, p in enumerate(per)))
print("\n== our emitted sign on pool genes, overall and where the donor consensus agrees/disagrees")
osg = np.sign(np.nan_to_num(d["ol"]))
a, n = acc(np.ones_like(pool), osg); print(f"  ours all: {a:.3f} (n {n})")
agree3 = (((nneg >= 2) & (npos == 0)) | ((npos >= 2) & (nneg == 0)))
a, n = acc(agree3 & (np.sign(cons) == osg), osg); print(f"  ours where >=2 donors agree WITH ours: {a:.3f} (n {n})")
a, n = acc(agree3 & (np.sign(cons) == -osg), osg); print(f"  ours where >=2 donors agree AGAINST ours: {a:.3f} (n {n})")
