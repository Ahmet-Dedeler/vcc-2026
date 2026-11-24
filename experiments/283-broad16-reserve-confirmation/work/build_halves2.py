"""Per working line, using the AUTHORS' per-cell assignment (processed/all_cells_table.csv: pass_qc, assigned_ko,
   assigned_grna): per-KO pseudobulk sums for FULL, seeded random cell halves A/B, and the two guides v1/v2, on the
   on-disk pseudobulk gene axis (duplicate gene names summed as the authors do).  Parse known-answer: FULL must equal
   the on-disk pseudobulk_sum_by_ko column exactly.   usage: build_halves2.py <line> <gefile> [seed]"""
import os, sys, time, json, numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from pathlib import Path
line, gefile = sys.argv[1:3]; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 268
t0 = time.time(); out = Path(os.environ.get("OUTDIR", "halves")); out.mkdir(exist_ok=True); EXCL = os.environ.get("EXCLUDE_ARM", "0") == "1"
d = np.load("pb_cache.npz", allow_pickle=True); X, genes, lines, targets = d["X"], d["genes"], d["lines"], d["targets"]
cols = np.where(lines == line)[0]; pb = pd.DataFrame(X[:, cols], index=genes, columns=targets[cols]); del X
tab = pd.read_csv("proc/all_cells_table.csv", low_memory=False); tab = tab[(tab.cell_line == line) & (tab.pass_qc == True)].copy()
if EXCL: tab = tab[(tab.arm_trunc == False) & (tab.arm_gain == False)].copy()
a = sc.read_10x_h5(gefile, gex_only=False)
bcs = np.array(a.obs_names); strip = "-" not in tab.barcode.iloc[0]
key = np.array([b.split("-")[0] for b in bcs]) if strip else bcs
pos = {b: i for i, b in enumerate(key)}; assert len(pos) == len(key), "barcode collision after stripping"
tab["row"] = tab.barcode.map(pos); n_missing = int(tab.row.isna().sum()); tab = tab.dropna(subset=["row"]); tab["row"] = tab.row.astype(int)
tab["ver"] = tab.assigned_grna.str.extract(r"_(v\d)_")[0]
# duplicate gene names summed (authors' anndata_groupby agg=sum)
vn = a.var_names.to_numpy(); uniq, inv = np.unique(vn, return_inverse=True)
Agg = sp.csr_matrix((np.ones(len(vn)), (np.arange(len(vn)), inv)), shape=(len(vn), len(uniq)))
M = (sp.csr_matrix(a.X, dtype=np.float32) @ Agg.astype(np.float32)).tocsr(); del a
gpos = {g: i for i, g in enumerate(uniq)}; idx = np.array([gpos.get(g, -1) for g in genes]); have = idx >= 0
lib = np.asarray(M.sum(1, dtype=np.float64)).ravel()
tlist = list(pb.columns); rng = np.random.default_rng(seed)
S = {k: np.full((len(tlist), len(genes)), np.nan, np.float32) for k in ("full", "A", "B", "v1", "v2")}
ncell = np.zeros((len(tlist), 3), int); libmean = np.zeros(len(tlist)); ka = []
groups = {k: v.row.to_numpy() for k, v in tab.groupby("assigned_ko")}; vers = {k: v.ver.to_numpy() for k, v in tab.groupby("assigned_ko")}
for j, t in enumerate(tlist):
    m = groups.get(t, np.array([], int)); tv = vers.get(t, np.array([]))
    if len(m) == 0: continue
    perm = rng.permutation(len(m)); A = m[perm[: len(m) // 2]]; B = m[perm[len(m) // 2:]]
    v1 = m[tv == "v1"]; v2 = m[tv == "v2"]
    for k, sel in (("full", m), ("A", A), ("B", B), ("v1", v1), ("v2", v2)):
        s = np.asarray(M[sel].sum(0, dtype=np.float64)).ravel() if len(sel) else np.zeros(len(uniq))
        row = np.full(len(genes), np.nan, np.float32); row[have] = s[idx[have]]; S[k][j] = row
    ncell[j] = (len(m), len(v1), len(v2)); libmean[j] = lib[m].mean()
    ref = pb[t].to_numpy(); mine = S["full"][j]; ok = ~np.isnan(ref) & ~np.isnan(mine)
    ka.append((t, bool(np.array_equal(ref[ok], mine[ok])), float(np.abs(ref[ok] - mine[ok]).sum() / max(ref[ok].sum(), 1)), int(ok.sum()), int((np.isnan(mine) & ~np.isnan(ref)).sum())))
cq = pd.read_csv("~/Code/virtual-cell-challenge-2026/data/external/depmap26q1/extracted/processed/cell_quality_covariates.csv").set_index("cell_line")
summary = {"line": line, "seed": seed, "exclude_arm": EXCL, "h5_cells": int(len(bcs)), "pass_qc_cells": int(len(tab) + n_missing), "unmatched_barcodes": n_missing, "kos_in_table": len(groups),
           "targets_with_cells": int((ncell[:, 0] > 0).sum()), "cells_per_ko_mean": float(ncell[:, 0][ncell[:, 0] > 0].mean()), "ref_cells_per_ko": float(cq.loc[line, "n_cells_per_ko"]),
           "cells_per_ko_min": int(ncell[:, 0][ncell[:, 0] > 0].min()), "guide_cells_min": int(ncell[ncell[:, 0] > 0, 1:].min()), "lib_mean": float(lib[tab.row].mean()), "ref_avg_umis": float(cq.loc[line, "avg_umis_per_ko"]),
           "ka_exact_targets": sum(k[1] for k in ka), "ka_targets": len(ka), "ka_median_rel_dev": float(np.median([k[2] for k in ka])), "ka_ref_genes_unmapped": int(max(k[4] for k in ka)), "genes_mapped": int(have.sum()), "seconds": round(time.time() - t0)}
print(json.dumps(summary), flush=True)
np.savez_compressed(out / f"{line}.npz", targets=np.array(tlist), genes=genes, ncell=ncell, libmean=libmean, **S)
json.dump({"summary": summary, "ka": ka}, open(out / f"{line}.json", "w"), indent=1)
