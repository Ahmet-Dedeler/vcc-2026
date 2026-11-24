"""515: does a low-rank (response-module) projection of a donor row give better per-gene SIGNS on the recipient's DE
pool than the raw row? H1 truth (h1-test, 85 targets), K562-GWPS rows (the only large single-lab screen that covers
H1's targets). The basis is learned on K562-GWPS perturbations that are NOT h1-test targets (no leakage of the rows being
projected). Also a Wiener-style shrink: row -> V diag(s2/(s2+n2)) V^T row with n2 the per-component noise floor
estimated from the tail singular values. Accuracy is measured on H1 pool genes (reference q < 0.05, finite, non-zero,
target excluded), overall and by pool band, only where both the raw and projected values are non-zero.
    uv run python experiments/515-lowrank-de-signs/lowrank515.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/513-reach-walk-anatomy"))
from walk513 import ALPHA, load  # noqa: E402

d = load(); G = d["genes"]; T = list(d["targets"])
B = ROOT / "data/donors/k562_gwps"
perts = [l.strip() for l in open(B / "perts.txt")]; genes = [l.strip() for l in open(B / "genes.txt")]
M = np.load(B / "delta.npy", mmap_mode="r")
gi = {g: i for i, g in enumerate(genes)}
common = [g for g in G if g in gi]
cols = np.array([gi[g] for g in common]); hcols = np.array([list(G).index(g) for g in common])
print(f"K562 genes on the H1 axis: {len(common)}")
train = [i for i, p in enumerate(perts) if p not in set(T)]
X = np.asarray(M[train][:, cols], dtype=np.float64)
X -= X.mean(0, keepdims=True) * 0  # no centring: the programme lives in the mean and we keep it
print(f"training rows {X.shape}")
# randomized SVD would do; the full thin SVD of ~9.7k x ~7k is fine in float64 (~0.5 GB)
U, s, Vt = np.linalg.svd(X, full_matrices=False)
print("singular values: top 5", np.round(s[:5], 2), " s[100]", round(s[100], 3), " s[1000]", round(s[1000], 3))
tail = np.median(s[len(s) // 2:]) ** 2          # per-component noise power proxy (bulk of the spectrum)
pool = np.zeros((len(T), len(common)), bool); truth = np.zeros((len(T), len(common)))
raw = np.full((len(T), len(common)), np.nan)
for k, t in enumerate(T):
    rl = d["rl"][k][hcols]
    pool[k] = (d["rq"][k][hcols] < ALPHA) & np.isfinite(rl) & (rl != 0) & (np.array(common) != t)
    truth[k] = np.sign(np.nan_to_num(rl))
    if t in perts:
        raw[k] = np.asarray(M[perts.index(t)][cols], dtype=np.float64)
has = np.isfinite(raw).all(1)
nr = pool.sum(1); band = np.where(nr < 20, 0, np.where(nr < 100, 1, 2))
def acc(P):
    out = []
    for b in (None, 0, 1, 2):
        m = pool & has[:, None] & (P != 0) & np.isfinite(P)
        if b is not None:
            m &= (band == b)[:, None]
        out.append(float((np.sign(P[m]) == truth[m]).mean()) if m.any() else np.nan)
    return out
rows = {"raw": raw}
R = np.nan_to_num(raw)
for k_ in (10, 30, 100, 300, 1000):
    V = Vt[:k_]
    rows[f"proj{k_}"] = np.where(has[:, None], (R @ V.T) @ V, np.nan)
w = s ** 2 / (s ** 2 + tail)
rows["wiener"] = np.where(has[:, None], ((R @ Vt.T) * w) @ Vt, np.nan)
rows["mean_prog"] = np.where(has[:, None], np.ones((len(T), 1)) * X.mean(0)[None, :], np.nan)   # the K562 programme alone
print(f"{'source':10s} {'all':>7s} {'<20':>7s} {'20-100':>7s} {'>=100':>7s}")
for name, P in rows.items():
    a = acc(P)
    print(f"{name:10s} " + " ".join(f"{x:7.3f}" for x in a))
