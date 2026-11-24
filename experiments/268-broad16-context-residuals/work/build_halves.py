"""Per working line: assign cells to targets from the guide-UMI table (rule: top guide >= min_top UMIs AND second guide
   <= max_second_frac * top), then write per-target pseudobulk sums for FULL, random cell halves A/B (seeded), and the two
   guides (v1/v2) on the pseudobulk gene axis, plus the parse known-answer vs the on-disk pseudobulk_sum_by_ko column.
usage: build_halves.py <line> <gefile> <guidefile> <min_top> <max_second_frac> [seed]"""
import sys, time, json, numpy as np, pandas as pd, h5py, scipy.sparse as sp
from pathlib import Path
line, gefile, gfile = sys.argv[1:4]; min_top = int(sys.argv[4]); max_sec = float(sys.argv[5]); seed = int(sys.argv[6]) if len(sys.argv) > 6 else 268
t0 = time.time(); out = Path("halves"); out.mkdir(exist_ok=True)
d = np.load("pb_cache.npz", allow_pickle=True); X, genes, lines, targets = d["X"], d["genes"], d["lines"], d["targets"]
cols = np.where(lines == line)[0]; pb = pd.DataFrame(X[:, cols], index=genes, columns=targets[cols]); del X
with h5py.File(gefile, "r") as h:
    bc = np.array([b.decode().split("-")[0] for b in h["matrix/barcodes"][:]]); vn = np.array([n.decode() for n in h["matrix/features/name"][:]])
    shape = tuple(h["matrix/shape"][:]); M = sp.csc_matrix((h["matrix/data"][:], h["matrix/indices"][:], h["matrix/indptr"][:]), shape=shape).tocsr()  # 10x stores genes x cells CSC -> cells x genes CSR
M = M.T.tocsr() if M.shape[0] != len(bc) else M
assert M.shape == (len(bc), len(vn)), (M.shape, len(bc), len(vn))
print(line, "h5", M.shape, "nnz", M.nnz, f"{time.time()-t0:.0f}s", flush=True)
hdr = pd.read_csv(gfile, nrows=0).columns; keepcols = [c for c in hdr if c.split("-")[0] in set(bc)]
g = pd.read_csv(gfile, index_col=0, usecols=[hdr[0]] + keepcols); g.columns = [c.split("-")[0] for c in g.columns]
gnames = np.array(g.index); gtarget = np.array([n.split("_")[0] for n in gnames]); gver = np.array([n.split("_")[1] for n in gnames])
pos = {b: i for i, b in enumerate(bc)}; rows = np.array([pos[c] for c in g.columns]); G = g.to_numpy(np.int64); del g
order = np.argsort(-G, axis=0); top = order[0]; second = order[1]; ar = np.arange(G.shape[1])
topc = G[top, ar]; secc = G[second, ar]
assigned = (topc >= min_top) & (secc <= max_sec * topc); tg = gtarget[top]; tv = gver[top]
dup = pd.Series(vn).duplicated(keep=False).to_numpy(); gpos = {n: i for i, n in enumerate(vn) if not dup[i]}
idx = np.array([gpos.get(gn, -1) for gn in genes]); have = idx >= 0
lib = np.asarray(M.sum(1, dtype=np.int64)).ravel()
tlist = list(pb.columns); rng = np.random.default_rng(seed)
S = {k: np.full((len(tlist), len(genes)), np.nan, np.float32) for k in ("full", "A", "B", "v1", "v2")}
ncell = np.zeros((len(tlist), 3), int); libmean = np.zeros(len(tlist)); ka = []
for j, t in enumerate(tlist):
    m = np.where(assigned & (tg == t))[0]
    if len(m) == 0: continue
    perm = rng.permutation(m); A = perm[: len(m) // 2]; B = perm[len(m) // 2:]
    v1 = m[tv[m] == "v1"]; v2 = m[tv[m] == "v2"]
    for k, sel in (("full", m), ("A", A), ("B", B), ("v1", v1), ("v2", v2)):
        s = np.asarray(M[rows[sel]].sum(0, dtype=np.int64)).ravel() if len(sel) else np.zeros(len(vn))
        row = np.full(len(genes), np.nan, np.float32); row[have] = s[idx[have]]; S[k][j] = row
    ncell[j] = (len(m), len(v1), len(v2)); libmean[j] = lib[rows[m]].mean()
    ref = pb[t].to_numpy(); mine = S["full"][j]; ok = ~np.isnan(ref) & ~np.isnan(mine)
    ka.append((t, bool(np.array_equal(ref[ok], mine[ok])), float(np.abs(ref[ok] - mine[ok]).sum() / max(ref[ok].sum(), 1)), float(np.corrcoef(ref[ok], mine[ok])[0, 1]), int(ok.sum())))
exact = sum(k[1] for k in ka); dev = float(np.median([k[2] for k in ka])); cor = float(np.median([k[3] for k in ka]))
cq = pd.read_csv("~/Code/virtual-cell-challenge-2026/data/external/depmap26q1/extracted/processed/cell_quality_covariates.csv").set_index("cell_line")
summary = {"line": line, "rule": {"min_top": min_top, "max_second_frac": max_sec, "seed": seed}, "ge_cells": int(len(bc)), "guide_cells": int(G.shape[1]), "assigned_cells": int(assigned.sum()),
           "targets_with_cells": int((ncell[:, 0] > 0).sum()), "cells_per_ko_mean": float(ncell[:, 0][ncell[:, 0] > 0].mean()), "ref_cells_per_ko": float(cq.loc[line, "n_cells_per_ko"]),
           "cells_per_ko_min": int(ncell[:, 0][ncell[:, 0] > 0].min()), "guide_min_cells": int(ncell[:, 1:].min()), "lib_mean": float(lib[rows[assigned]].mean()), "ref_avg_umis": float(cq.loc[line, "avg_umis_per_ko"]),
           "ka_exact_targets": exact, "ka_median_rel_dev": dev, "ka_median_corr": cor, "ka_targets": len(ka), "genes_mapped": int(have.sum()), "seconds": round(time.time() - t0)}
print(json.dumps(summary), flush=True)
np.savez_compressed(out / f"{line}.npz", targets=np.array(tlist), genes=genes, ncell=ncell, libmean=libmean, **S)
json.dump({"summary": summary, "ka": ka}, open(out / f"{line}.json", "w"), indent=1)
