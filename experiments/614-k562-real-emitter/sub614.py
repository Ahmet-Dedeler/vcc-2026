#!/usr/bin/env python3
"""614 amendment A1: the K562 truth plus an H1-shaped SUB-GATE RESERVOIR.
The production emitter balances every emitted cell's exact library with sub-gate genes (< 5 CPM, outside the DE
universe). The GWPS release keeps only 8,248 expressed genes, so the K562 truth has none, and rc_hybrid_k cannot meet the
per-cell totals (every K562 run crashed). This appends 7,301 synthetic sub-gate genes, mirroring H1's reservoir:
  - its count (7,301 = h1-test's non-kept genes);
  - its share of each cell's library (s0 = 0.0058 = H1's control sub-gate mass);
  - its per-gene profile (gene j's share = H1 sub-gate gene j's share of H1's sub-gate mass), so each gene sits at H1's
    sub-5-CPM level.
Per cell i: T_i ~ Poisson(L_i · s0 / (1 − s0)), spread multinomially over the 7,301 genes. Control and perturbed cells
are drawn identically (no effect), seed 614. DE members never read these genes. They change every cell's CPM
normaliser by the same ~0.58 %, as H1's real sub-gate genes do.
    uv run python experiments/614-k562-real-emitter/sub614.py  -> work/k562-gwps-c250-sub.h5"""
import json
import os
from pathlib import Path
import h5py
import numpy as np
import scipy.sparse as sp
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
r = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
kh = np.asarray(r["keep"], bool)
cs = np.asarray(r["sums"], float)[-1]
w = cs[~kh] / cs[~kh].sum()
SCALE = float(os.environ.get("SUB614_SCALE", "1"))   # A1-R: reservoir x0.5 / x2 sensitivity (default 1 = A1 itself)
s0 = float(cs[~kh].sum() / cs.sum()) * SCALE
NF = int((~kh).sum())
src, dst = HERE / "work/k562-gwps-c250.h5", HERE / ("work/k562-gwps-c250-sub.h5" if SCALE == 1 else f"work/k562-gwps-c250-sub{SCALE:g}.h5")
with h5py.File(src, "r") as h:
    lab = h["obs/target_gene"][:]
    genes = h["var/_index"][:].astype(str)
    shape = tuple(int(v) for v in h["X"].attrs["shape"])
    X = sp.csr_matrix((h["X/data"][:], h["X/indices"][:], h["X/indptr"][:]), shape=shape)
L = np.asarray(X.sum(1)).ravel()
rng = np.random.default_rng(614)
T = rng.poisson(L * s0 / (1 - s0))
rows, cols, vals = [], [], []
for i in range(shape[0]):
    if T[i] == 0:
        continue
    c = rng.multinomial(int(T[i]), w)
    nz = np.flatnonzero(c)
    rows.append(np.full(nz.size, i)); cols.append(nz); vals.append(c[nz])
S = sp.csr_matrix((np.concatenate(vals).astype(np.float32), (np.concatenate(rows), np.concatenate(cols))), shape=(shape[0], NF))
Y = sp.hstack([X, S], format="csr")
Y.sort_indices()
names = np.concatenate([genes, np.array([f"SUBGATE_{j:05d}" for j in range(NF)])])
with h5py.File(dst, "w") as o:
    o.create_dataset("obs/target_gene", data=lab)
    o.create_dataset("var/_index", data=names.astype("S"))
    x = o.create_group("X")
    x.attrs["shape"] = Y.shape
    x.create_dataset("data", data=Y.data.astype(np.float32))
    x.create_dataset("indices", data=Y.indices.astype(np.int32))
    x.create_dataset("indptr", data=Y.indptr.astype(np.int64))
ctrl = lab.astype(str) == "non-targeting"
L2 = np.asarray(Y.sum(1)).ravel()
mc = 1e6 * np.asarray(S[ctrl].multiply(1.0 / L2[ctrl][:, None]).mean(0)).ravel()
info = {"n_synthetic": NF, "share_s0": s0, "median_control_library_before": float(np.median(L[ctrl])),
        "median_control_library_after": float(np.median(L2[ctrl])), "synthetic_mean_cpm_max": float(mc.max()),
        "synthetic_genes_over_5cpm": int((mc > 5).sum()), "nnz_added": int(S.nnz)}
(HERE / ("work/sub614.json" if SCALE == 1 else f"work/sub614-x{SCALE:g}.json")).write_text(json.dumps(info, indent=1))
print(json.dumps(info))
