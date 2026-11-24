"""434e: the top of OUR per-target ranking. Per target, genes ranked by our -log10 q (then |lfc|). On h1-test's gates:
sign accuracy and gate-mass share by rank bucket, split by target gate size; then the rule 'deliver magnitude M on our
top-k genes (our sign), plain MDE elsewhere' for a grid of k and M. Analytic nmae only."""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
z = np.load(HERE / "work/allcells434.npz", allow_pickle=True)
L, Q, T, D = z["L"].astype(float), z["Q"].astype(float), z["T"].astype(float), z["D"].astype(float)
donors = [str(x) for x in z["donors"]]
r = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
keep = r["keep"].astype(bool); ref_lin = np.asarray(r["ref_lin"], float)
ref_lin = ref_lin[keep] if ref_lin.shape[0] == keep.shape[0] else ref_lin
cnt = ref_lin * np.median(np.asarray(r["libraries"], float)) / 1e6
S = np.sign(L)
rank = np.zeros_like(L, dtype=int)
for k in range(L.shape[0]):
    key = np.lexsort((-np.abs(L[k]), np.log10(Q[k] + 1e-300)))       # q ascending, then |lfc| descending
    rank[k, key] = np.arange(L.shape[1])
MDE = np.where(cnt >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(cnt, 1e-3)), 0.1), 0.0)[None, :]
base_pred = L + S * MDE
tg = [k for k in range(L.shape[0]) if np.isfinite(T[k]).any()]
gsize = {k: int(np.isfinite(T[k]).sum()) for k in tg}

def nmae(P, ks=tg):
    return float(np.mean([np.abs(P[k, np.isfinite(T[k])] - T[k, np.isfinite(T[k])]).sum() / np.abs(T[k, np.isfinite(T[k])]).sum() for k in ks]))

buckets = [(0, 10), (10, 50), (50, 200), (200, 1000), (1000, 3000), (3000, 99999)]
for nm, ks in (("all", tg), ("gate<100", [k for k in tg if gsize[k] < 100]), ("gate>=100", [k for k in tg if gsize[k] >= 100])):
    print(f"\n{nm} ({len(ks)} targets): rank bucket | gate cells | sign acc | |l|-w acc | share of gate |l| | median |l| | median |ours|")
    tot = sum(np.abs(T[k][np.isfinite(T[k])]).sum() for k in ks)
    for lo, hi in buckets:
        cells = [(k, g) for k in ks for g in np.flatnonzero(np.isfinite(T[k]) & (rank[k] >= lo) & (rank[k] < hi))]
        if not cells:
            continue
        kk, gg = np.array(cells).T
        l = T[kk, gg]; s = S[kk, gg]
        print(f"   [{lo:5d},{hi:5d}) | {len(cells):6d} | {np.mean(s == np.sign(l)):.3f} | {np.sum((s == np.sign(l)) * np.abs(l)) / np.abs(l).sum():.3f} | "
              f"{np.abs(l).sum() / tot:.3f} | {np.median(np.abs(l)):.3f} | {np.median(np.abs(L[kk, gg])):.3f}")
b0 = nmae(L); bm = nmae(base_pred)
print(f"\nbase {b0:.5f}  MDE {bm - b0:+.5f}")
for k_ in (5, 10, 20, 50, 100):
    row = []
    for M in (0.15, 0.25, 0.4, 0.6):
        top = rank < k_
        P = np.where(top, S * np.maximum(np.abs(L), M), base_pred)
        row.append(nmae(P) - b0)
    print(f"  top-{k_:3d}: M 0.15 {row[0]:+.5f}  0.25 {row[1]:+.5f}  0.4 {row[2]:+.5f}  0.6 {row[3]:+.5f}")
# panel-donor-agreement gated version: top-k AND >=2 of the panel donors agree with our sign
PAN = [donors.index(d) for d in ("k562_gwps", "orion_hct116", "orion_hek293t")]
agree = np.nansum(np.sign(D[PAN]) * S[None], 0)
for k_ in (10, 20, 50):
    row = []
    for M in (0.25, 0.4, 0.6):
        top = (rank < k_) & (agree >= 2)
        P = np.where(top, S * np.maximum(np.abs(L), M), base_pred)
        row.append(nmae(P) - b0)
    print(f"  top-{k_:3d} & 2 panel donors agree: M 0.25 {row[0]:+.5f}  0.4 {row[1]:+.5f}  0.6 {row[2]:+.5f}")
