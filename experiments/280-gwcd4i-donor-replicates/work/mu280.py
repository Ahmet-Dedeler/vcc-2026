"""280 -- Q1's headline in the AUTHORS' OWN estimator: `GWCD4i.DE_stats.by_donors.h5mu` holds a DESeq2 fit
computed inside each of the 6 donor pairs.  The three disjoint pairings give exactly the §3.2 statistic with
DESeq2 log2 fold changes instead of a log1p(CPM) pseudobulk difference.

Read over HTTP range requests -- the `layers/log_fc` datasets are CONTIGUOUS and uncompressed, so only the
panel rows are fetched (nothing is downloaded to disk).
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import h5py

sys.path.insert(0, ".")
from lib280 import boot_median, boot_median_diff, cosines  # noqa: E402
from remote_h5 import RemoteFile  # noqa: E402

ROOT = "~/Code/virtual-cell-challenge-2026"
B = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
P = lambda *a: print(*a, flush=True)  # noqa: E731
D = {"CE0008162": "D1", "CE0010866": "D2", "CE0008678": "D3", "CE0006864": "D4"}
PAIRS = {  # the three DISJOINT pairings, as modality names
    "D1D2|D3D4": ("CE0008162_CE0010866", "CE0006864_CE0008678"),
    "D1D3|D2D4": ("CE0008162_CE0008678", "CE0006864_CE0010866"),
    "D1D4|D2D3": ("CE0006864_CE0008162", "CE0008678_CE0010866"),
}
panel = sorted(set(l.split(",")[0] for l in
                   open(f"{ROOT}/data/challenge/pert_counts.csv").read().strip().split("\n")[1:]))
panel_set = set(panel)
t0 = time.time()
rf = RemoteFile(f"{B}/GWCD4i.DE_stats.by_donors.h5mu", block=1 << 20, max_blocks=1024)
f = h5py.File(rf, "r")
mods = sorted(f["mod"].keys())
P(f"{len(mods)} modalities: {', '.join(D[m.split('_')[0]] + '+' + D[m.split('_')[1]] for m in mods)}")

store: dict = {}
for m in mods:
    g = f[f"mod/{m}"]
    def cat(n):
        gg = g[f"obs/{n}"]
        return np.asarray(gg["categories"]).astype(str)[np.asarray(gg["codes"])]
    tgt, cc = cat("target_contrast_gene_name"), cat("culture_condition")
    ncell = np.asarray(g["obs/n_cells_target"])
    genes = np.asarray(g["var/gene_name"]).astype(str)
    sel = np.where(np.isin(tgt, panel))[0]
    lf = g["layers/log_fc"]
    X = np.empty((len(sel), lf.shape[1]))
    for j, i in enumerate(sel):
        X[j] = lf[int(i)]
    store[m] = dict(target=tgt[sel], cond=cc[sel], ncell=ncell[sel], genes=genes, lfc=X)
    P(f"  {D[m.split('_')[0]]}+{D[m.split('_')[1]]}: {lf.shape[0]} obs, {len(sel)} panel rows fetched"
      f"   [{time.time()-t0:.0f}s, {rf.n_requests} req, {rf.n_bytes/1e9:.2f} GB]")

# ---- gene axis: the authors' 10,273 genes -> the official axis, panel genes dropped
official = [l.strip() for l in open(f"{ROOT}/data/challenge/gene_names.csv").read().strip().split("\n")[1:]]
offpos = {gname: i for i, gname in enumerate(official)}
genes = store[mods[0]]["genes"]
assert all((store[m]["genes"] == genes).all() for m in mods), "modalities differ in var order"
keep = np.asarray([(gname in offpos) and (gname not in panel_set) for gname in genes])
P(f"gene axis: {len(genes)} authors' genes -> {int(keep.sum())} scored (on the official axis, non-panel)")

rows, per = [], {}
for cond in ("Rest", "Stim8hr", "Stim48hr"):
    store_c = {}
    for m in mods:
        s = store[m]
        k = s["cond"] == cond
        store_c[m] = pd.DataFrame(dict(target=s["target"][k], ncell=s["ncell"][k],
                                       row=np.where(k)[0]))
    acc = {}
    for lab, (a, b) in PAIRS.items():
        common = sorted(set(store_c[a].target) & set(store_c[b].target))
        ia = store_c[a].set_index("target").loc[common, "row"].values
        ib = store_c[b].set_index("target").loc[common, "row"].values
        A = store[a]["lfc"][ia][:, keep]
        Bm = store[b]["lfc"][ib][:, keep]
        A = A - A.mean(0, keepdims=True)          # the same panel residual centring as §3.2
        Bm = Bm - Bm.mean(0, keepdims=True)
        cs = cosines(A, Bm)
        acc[lab] = pd.Series(cs, index=common)
        m_, lo, hi, n = boot_median(cs)
        P(f"  {cond:<9} {lab}  n={n:3d}  median {m_:+.4f} [{lo:+.4f},{hi:+.4f}]  mean {np.nanmean(cs):+.4f}"
          f"  <=0 {np.nanmean(cs<=0):.3f}")
        rows.append(dict(cond=cond, pairing=lab, n=n, median=m_, lo=lo, hi=hi,
                         mean=float(np.nanmean(cs)), frac_le0=float(np.nanmean(cs <= 0))))
    allt = sorted(set().union(*[set(v.index) for v in acc.values()]))
    pooled = pd.DataFrame({k: v for k, v in acc.items()}).reindex(allt).mean(1)
    m_, lo, hi, n = boot_median(pooled.values)
    P(f"  {cond:<9} POOLED (mean over the 3 pairings): n={n} median {m_:+.4f} [{lo:+.4f},{hi:+.4f}]")
    rows.append(dict(cond=cond, pairing="POOLED", n=n, median=m_, lo=lo, hi=hi,
                     mean=float(np.nanmean(pooled.values)), frac_le0=float(np.nanmean(pooled.values <= 0))))
    per[cond] = pooled

pd.DataFrame(rows).to_csv("mu280_donor_pairs.csv", index=False)

# ---- agreement with my pseudobulk statistic, and with 276
P("\n=== agreement with the log1p(CPM) pseudobulk statistic of §2.1, and with 276 ===")
c276 = pd.read_csv(f"{ROOT}/experiments/276-conservation-map/work/conservation276.csv").set_index("target")
out = {}
for cond in ("Rest", "Stim8hr", "Stim48hr"):
    z = np.load(f"q1_pertarget_{cond}_keepDE.npz", allow_pickle=True)
    mine = pd.Series(z["cos"], index=[str(x) for x in z["targets"]])
    j = pd.concat([per[cond].rename("deseq2"), mine.rename("pseudobulk")], axis=1).dropna()
    rho, p = spearmanr(j.deseq2, j.pseudobulk)
    v = c276.reindex(j.index)["C_raw"].values.astype(float)
    mk = np.isfinite(v)
    rho2, p2 = spearmanr(v[mk], j.deseq2.values[mk])
    P(f"  {cond:<9} Spearman(DESeq2 donor-pair cosine, my pseudobulk cosine) = {rho:+.3f} (p={p:.2g}, n={len(j)})"
      f"   |  Spearman(DESeq2, C_raw) = {rho2:+.3f} (p={p2:.2g}, n={mk.sum()})")
    out[cond] = dict(rho_vs_pseudobulk=float(rho), p=float(p), n=int(len(j)),
                     rho_vs_C_raw=float(rho2), p_C_raw=float(p2))
    per[cond].to_csv(f"mu280_pertarget_{cond}.csv", header=["cos"])
json.dump(out, open("mu280.json", "w"), indent=1)
P(f"\n[{time.time()-t0:.0f}s, {rf.n_requests} range requests, {rf.n_bytes/1e9:.2f} GB fetched, 0 B written to disk]")
