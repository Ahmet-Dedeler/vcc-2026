"""Gene-axis mapping CD4 pseudobulk (18,129) -> official challenge axis (18,533), and the
authors' own cross-donor correlations from the LOCAL DE_stats.h5ad (KA/cross-check material)."""
import numpy as np, h5py, json, collections
ROOT = "~/Code/virtual-cell-challenge-2026"
z = np.load("obs280.npz", allow_pickle=True)
cd4_genes = z["genes"].astype(str)
official = [l.strip() for l in open(f"{ROOT}/data/challenge/gene_names.csv").read().strip().split("\n")[1:]]
panel = [l.split(",")[0] for l in open(f"{ROOT}/data/challenge/pert_counts.csv").read().strip().split("\n")[1:]]
print("official", len(official), "unique", len(set(official)), "| cd4", len(cd4_genes), "unique", len(set(cd4_genes)))
dup = [g for g, c in collections.Counter(cd4_genes).items() if c > 1]
print("cd4 duplicated symbols:", len(dup), dup[:8])
pos = {}
for i, g in enumerate(cd4_genes):
    pos.setdefault(g, []).append(i)
cols_off, cols_cd4 = [], []
for j, g in enumerate(official):
    p = pos.get(g)
    if p is not None and len(p) == 1:
        cols_off.append(j); cols_cd4.append(p[0])
cols_off = np.asarray(cols_off); cols_cd4 = np.asarray(cols_cd4)
panel_set = set(panel)
is_panel = np.asarray([official[j] in panel_set for j in cols_off])
print(f"mapped genes: {len(cols_off)} of official {len(official)}; of which panel genes {is_panel.sum()}")
print(f"SCORED columns (mapped, non-panel): {(~is_panel).sum()}")
np.savez("axis280.npz", cols_off=cols_off, cols_cd4=cols_cd4, is_panel=is_panel,
         official=np.asarray(official), panel=np.asarray(panel))

# ---- the authors' own donor correlations, from the local DE table
f = h5py.File(f"{ROOT}/experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad", "r")
def cat(n):
    g = f["obs/" + n]; c = np.asarray(g["categories"]).astype(str); return c[np.asarray(g["codes"])]
tg, cc = cat("target_contrast_gene_name"), cat("culture_condition")
cols = {k: np.asarray(f["obs/" + k]) for k in ("n_cells_target", "n_guides", "donor_correlation_all_mean",
        "donor_correlation_all_min", "donor_correlation_hits_mean", "donor_correlation_hits_min",
        "guide_correlation_all", "guide_correlation_signif", "n_total_de_genes", "ontarget_effect_size")}
m = np.isin(tg, panel)
import pandas as pd
de = pd.DataFrame(dict(target=tg[m], cond=cc[m], **{k: v[m] for k, v in cols.items()}))
de.to_csv("de_stats_panel280.csv", index=False)
print("\n== authors' donor_correlation_all_mean on the 300-panel rows (per condition) ==")
print(de.groupby("cond").donor_correlation_all_mean.describe()[["count","mean","25%","50%","75%"]].to_string())
print("\n== authors' guide_correlation_all (same rows) ==")
print(de.groupby("cond").guide_correlation_all.describe()[["count","mean","25%","50%"]].to_string())
print("\n== n_cells_target ==")
print(de.groupby("cond").n_cells_target.describe()[["count","mean","25%","50%","75%","max"]].to_string())
