"""Read pseudobulk_merged obs + indptr + var over range requests; report the row layout and the
exact byte cost of the panel-target + NTC subset."""
import time, json
import numpy as np, h5py
from remote_h5 import RemoteFile
import pandas as pd

B = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
ROOT = "~/Code/virtual-cell-challenge-2026"
t0 = time.time()
rf = RemoteFile(f"{B}/GWCD4i.pseudobulk_merged.h5ad", block=1 << 22, max_blocks=64)
f = h5py.File(rf, "r")

def cat(name):
    g = f[f"obs/{name}"]
    c = np.asarray(g["categories"]).astype(str)
    return c[np.asarray(g["codes"])], c

donor, dcats = cat("donor_id")
cond, ccats = cat("culture_condition")
gtype, gtcats = cat("guide_type")
run, rcats = cat("10xrun_id")
pert, pcats = cat("perturbed_gene_name")
gid, gidcats = cat("guide_id")
ncells = np.asarray(f["obs/n_cells"])
tot = np.asarray(f["obs/total_counts"])
keep_de = np.asarray(f["obs/keep_for_DE"])
keep_min = np.asarray(f["obs/keep_min_cells"])
keep_eff = np.asarray(f["obs/keep_effective_guides"])
keep_tc = np.asarray(f["obs/keep_total_counts"])
keep_tg = np.asarray(f["obs/keep_test_genes"])
genes = np.asarray(f["var/gene_name"]).astype(str)
gene_ids = np.asarray(f["var/gene_ids"]).astype(str)
print("donor cats", list(dcats), "| cond", list(ccats), "| guide_type", list(gtcats), "| run", list(rcats))
print("n_obs", len(donor), "n_var", len(genes), f"[{time.time()-t0:.1f}s {rf.n_requests} req {rf.n_bytes/1e6:.0f} MB]")

df = pd.DataFrame(dict(donor=donor, cond=cond, gtype=gtype, run=run, pert=pert, guide=gid,
                       n_cells=ncells, total=tot, keep_for_DE=keep_de, keep_min=keep_min,
                       keep_eff=keep_eff, keep_tc=keep_tc, keep_tg=keep_tg))
df["row"] = np.arange(len(df))
print("\n=== rows per donor x condition x guide_type ===")
print(df.groupby(["cond", "donor", "gtype"]).size().unstack(fill_value=0))
print("\n=== donor x run per condition ===")
print(df.groupby(["cond", "donor", "run"]).size().unstack(fill_value=0))

panel = [l.split(",")[0] for l in open(f"{ROOT}/data/challenge/pert_counts.csv").read().strip().split("\n")[1:]]
panel = sorted(set(panel))
print("\npanel targets", len(panel), "| present in pseudobulk pert names:", len(set(panel) & set(pcats)))
sel = df[(df.pert.isin(panel)) | (df.gtype != "targeting")]
print("selected rows (panel targeting + all non-targeting):", len(sel))
print(sel.groupby(["cond", "donor", "gtype"]).size().unstack(fill_value=0))

indptr = np.asarray(f["X/indptr"])
nnz = indptr[1:] - indptr[:-1]
rows = np.sort(sel.row.values)
b_data = int(nnz[rows].sum()) * 8
b_idx = int(nnz[rows].sum()) * 8
print(f"\nnnz total {indptr[-1]:,}; selected nnz {nnz[rows].sum():,} -> data+indices = {(b_data+b_idx)/1e9:.2f} GB")
runs = np.split(rows, np.where(np.diff(rows) != 1)[0] + 1)
print(f"contiguous row runs: {len(runs)} (max len {max(len(r) for r in runs)})")
np.savez_compressed("obs280.npz", donor=donor, cond=cond, gtype=gtype, run=run, pert=pert, guide=gid,
                    n_cells=ncells, total=tot, keep_for_DE=keep_de, keep_min=keep_min, keep_eff=keep_eff,
                    keep_tc=keep_tc, keep_tg=keep_tg, genes=genes, gene_ids=gene_ids, indptr=indptr)
print(f"\n[{time.time()-t0:.1f}s total, {rf.n_requests} range requests, {rf.n_bytes/1e6:.0f} MB fetched]")
