"""Gene-level comparison of two variants' nmae contributions over all scored targets (stats npz)."""
import sys, numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN
z = np.load(RN.ROOT / "experiments/323-current-family-validation/work/reference.npz")
genes_k = z["genes"][z["keep"]]; targets = z["targets"]; RQ = z["q"]; RL = z["lfc"]
tag = sys.argv[3] if len(sys.argv) > 3 else "s0"
B = np.load(HERE / f"work/h1/stats-base_{tag}.npz")["lfc"]
V1 = np.load(HERE / f"work/h1/stats-{sys.argv[1]}_{tag}.npz")["lfc"]; V2 = np.load(HERE / f"work/h1/stats-{sys.argv[2]}_{tag}.npz")["lfc"]
cat = {"both": [0, 0, 0], "only1": [0, 0, 0], "only2": [0, 0, 0]}
recs = []
nt = 0
for i in range(len(targets)):
    S = (RQ[i] < .05) & np.isfinite(RL[i]) & (RL[i] != 0) & (genes_k != targets[i])
    if S.sum() < 10: continue
    nt += 1
    r = RL[i]; den = np.abs(r[S]).sum(); b = np.nan_to_num(B[i]); a1 = np.nan_to_num(V1[i]); a2 = np.nan_to_num(V2[i])
    c1 = np.where(S, (np.abs(a1 - r) - np.abs(b - r)) / den, 0); c2 = np.where(S, (np.abs(a2 - r) - np.abs(b - r)) / den, 0)
    m1 = S & (a1 != b); m2 = S & (a2 != b)
    for k, m in (("both", m1 & m2), ("only1", m1 & ~m2), ("only2", ~m1 & m2)):
        cat[k][0] += m.sum(); cat[k][1] += c1[m].sum(); cat[k][2] += c2[m].sum()
    for j in np.flatnonzero(m1 | m2):
        recs.append((c1[j] - c2[j], targets[i], genes_k[j], r[j], b[j], a1[j], a2[j]))
for k, v in cat.items():
    print(f"{k}: genes {v[0]}, sum contrib {sys.argv[1]} {v[1] / nt:+.5f}, {sys.argv[2]} {v[2] / nt:+.5f} (per-target mean)")
recs.sort()
print("largest advantages of", sys.argv[1])
for rr in recs[:10]: print("  %s %s ref %+.2f base %+.3f v1 %+.3f v2 %+.3f d %+.5f" % (rr[1], rr[2], rr[3], rr[4], rr[5], rr[6], rr[0] / nt))
print("largest advantages of", sys.argv[2])
for rr in recs[-8:]: print("  %s %s ref %+.2f base %+.3f v1 %+.3f v2 %+.3f d %+.5f" % (rr[1], rr[2], rr[3], rr[4], rr[5], rr[6], rr[0] / nt))
