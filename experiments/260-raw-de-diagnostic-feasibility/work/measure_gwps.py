"""K562-GWPS availability numbers from obs/var only, plus ONE X row-block sample.
Never loads X. Peak RSS well under 1 GB."""
import json, time, h5py, numpy as np, pandas as pd
from pathlib import Path
P = "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
OUT = Path("experiments/260-raw-de-diagnostic-feasibility/work")
panel = pd.read_csv("data/challenge/pert_counts.csv")["target_gene"].astype(str).tolist()
axis = pd.read_csv("data/challenge/gene_names.csv")
axis_genes = axis.iloc[:, 0].astype(str).tolist() if axis.shape[1] == 1 else axis["gene_name"].astype(str).tolist()
print("panel", len(panel), "axis genes", len(axis_genes), "axis cols", list(axis.columns))

def cat(f, key):
    g = f[key]; cats = g["categories"][:].astype(str); codes = g["codes"][:]
    return cats, codes

res = {}
with h5py.File(P, "r") as f:
    o = f["obs"]
    cats, codes = cat(o, "perturbation")
    gcats, gcodes = cat(o, "gene")
    print("perturbation cats", len(cats), "gene cats", len(gcats), "same labels:",
          set(cats) == set(gcats))
    ncounts = o["ncounts"][:]; umi = o["UMI_count"][:]; batch = o["batch"][:]
    nperts = o["nperts"][:]; ngenes = o["ngenes"][:]
    print("ncounts==UMI_count allclose:", np.allclose(ncounts, umi), "nperts uniq", np.unique(nperts)[:5])
    counts = np.bincount(codes, minlength=len(cats))
    ctrl_i = list(cats).index("control")
    print("control cells", counts[ctrl_i])
    print("median UMI all cells", np.median(ncounts), "control", np.median(ncounts[codes == ctrl_i]),
          "pert", np.median(ncounts[codes != ctrl_i]))
    q = np.percentile(ncounts, [5, 25, 50, 75, 95])
    print("UMI pctiles 5/25/50/75/95", q)
    res["median_umi_all"] = float(np.median(ncounts)); res["median_umi_ctrl"] = float(np.median(ncounts[codes == ctrl_i]))
    res["median_umi_pert"] = float(np.median(ncounts[codes != ctrl_i])); res["umi_pctiles_5_25_50_75_95"] = q.tolist()
    res["n_cells"] = int(len(codes)); res["n_control_cells"] = int(counts[ctrl_i]); res["n_groups_excl_control"] = int(len(cats) - 1)
    res["median_genes_per_cell"] = float(np.median(ngenes))
    # batches
    ub, bc = np.unique(batch, return_counts=True)
    print("n batches", len(ub), "cells/batch min/median/max", bc.min(), np.median(bc), bc.max())
    bmed = np.array([np.median(ncounts[batch == b]) for b in ub])
    print("per-batch median UMI min/median/max", bmed.min(), np.median(bmed), bmed.max())
    ctrl_frac = np.array([np.mean(codes[batch == b] == ctrl_i) for b in ub])
    print("per-batch control fraction min/median/max", ctrl_frac.min(), np.median(ctrl_frac), ctrl_frac.max())
    res["n_batches"] = int(len(ub)); res["batch_median_umi_min_med_max"] = [float(bmed.min()), float(np.median(bmed)), float(bmed.max())]
    res["batch_ctrl_frac_min_med_max"] = [float(ctrl_frac.min()), float(np.median(ctrl_frac)), float(ctrl_frac.max())]
    # all-group cell count distribution
    pc = np.delete(counts, ctrl_i)
    print("cells/group: min/median/max", pc.min(), np.median(pc), pc.max(),
          ">=150:", (pc >= 150).sum(), ">=400:", (pc >= 400).sum(), "of", len(pc))
    res["all_groups_ge150"] = int((pc >= 150).sum()); res["all_groups_ge400"] = int((pc >= 400).sum())
    # panel overlap
    idx = {c: i for i, c in enumerate(cats)}
    inK = [t for t in panel if t in idx]
    print("panel targets in GWPS:", len(inK), "/", len(panel))
    pcnt = {t: int(counts[idx[t]]) for t in inK}
    arr = np.array(list(pcnt.values()))
    print("panel per-target cells: min/25/50/75/max", arr.min(), *np.percentile(arr, [25, 50, 75]).tolist(), arr.max())
    for th in (100, 150, 200, 250, 300, 400):
        print(f"  panel targets with >= {th} cells: {(arr >= th).sum()}")
    res["panel_in_gwps"] = len(inK); res["panel_cells_quartiles"] = [float(x) for x in np.percentile(arr, [0, 25, 50, 75, 100])]
    res["panel_ge"] = {str(th): int((arr >= th).sum()) for th in (100, 150, 200, 250, 300, 400)}
    # per-target median UMI for panel targets vs control
    pmed = {t: float(np.median(ncounts[codes == idx[t]])) for t in inK}
    ratio = np.array([pmed[t] for t in inK]) / res["median_umi_ctrl"]
    print("panel per-target depth / control: median", np.median(ratio), "IQR", np.percentile(ratio, [25, 75]))
    res["panel_depth_ratio_median_iqr"] = [float(np.median(ratio)), *[float(x) for x in np.percentile(ratio, [25, 75])]]
    pd.DataFrame({"target": inK, "cells": [pcnt[t] for t in inK], "median_umi": [pmed[t] for t in inK]}).sort_values("cells", ascending=False).to_csv(OUT / "panel_targets_in_gwps.csv", index=False)
    missing = [t for t in panel if t not in idx]
    res["panel_missing"] = missing
    print("panel targets missing from GWPS:", len(missing), missing[:40])
    # gene axis
    vg = f["var"]["gene_name"][:].astype(str)
    inter = sorted(set(vg) & set(axis_genes))
    print("GWPS genes", len(vg), "in official axis:", len(inter), "; panel targets among GWPS genes:", len(set(vg) & set(panel)))
    res["n_genes_gwps"] = int(len(vg)); res["n_genes_in_official_axis"] = len(inter); res["panel_targets_in_gwps_gene_axis"] = int(len(set(vg) & set(panel)))
    # ONE X block: integrality + rate
    X = f["X"]; print("X", X.shape, X.dtype, "chunks", X.chunks, "compression", X.compression, X.compression_opts)
    t0 = time.time(); blk = X[0:3886, :]; dt = time.time() - t0
    rs = blk.sum(1)
    print(f"block read {blk.shape} in {dt:.1f}s -> full pass est {dt * (X.shape[0] / 3886) / 60:.1f} min")
    print("integral:", bool(np.all(np.mod(blk, 1) == 0)), "min", blk.min(), "max", blk.max(),
          "rowsum==ncounts:", np.allclose(rs, ncounts[0:3886]), "rowsum==UMI_count:", np.allclose(rs, umi[0:3886]),
          "rowsum median", np.median(rs), "ncounts median (block)", np.median(ncounts[0:3886]))
    res["x_integral_block"] = bool(np.all(np.mod(blk, 1) == 0)); res["block_read_s"] = dt
    res["full_pass_est_min"] = dt * (X.shape[0] / 3886) / 60
    res["rowsum_eq_ncounts"] = bool(np.allclose(rs, ncounts[0:3886])); res["rowsum_eq_umi_count"] = bool(np.allclose(rs, umi[0:3886]))
    res["block_rowsum_median"] = float(np.median(rs)); res["block_ncounts_median"] = float(np.median(ncounts[0:3886]))
json.dump(res, open(OUT / "gwps_numbers.json", "w"), indent=1)
# essential ∩ panel
ess = set(open("data/donors/k562_essential/perts.txt").read().split())
print("k562_essential perts", len(ess), "∩ panel", len(ess & set(panel)))
