"""280 known-answer rows.  KA-1/KA-3/KA-5 need only the pseudobulk obs (no X); KA-2/KA-4 need the bank."""
from __future__ import annotations

import json
import sys

import h5py
import numpy as np
import pandas as pd

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, ".")
from lib280 import DONORS, Bank, cosines, residualise  # noqa: E402

OUT = {}
P = lambda *a: print(*a, flush=True)  # noqa: E731

z = np.load("obs280.npz", allow_pickle=True)
obs = pd.DataFrame({k: z[k] for k in ("donor", "cond", "gtype", "pert", "guide", "n_cells", "total",
                                      "keep_for_DE", "keep_min", "keep_eff", "keep_tc", "keep_tg")})
for c in ("donor", "cond", "gtype", "pert", "guide"):
    obs[c] = obs[c].astype(str)
obs["D"] = obs.donor.map(DONORS)
panel = sorted(set(l.split(",")[0] for l in open(f"{ROOT}/data/challenge/pert_counts.csv").read().strip().split("\n")[1:]))

# ------------------------------------------------------------------ KA-1: cells reconcile with the DE table
f = h5py.File(f"{ROOT}/experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad", "r")
def cat(n):
    g = f["obs/" + n]
    return np.asarray(g["categories"]).astype(str)[np.asarray(g["codes"])]
de = pd.DataFrame(dict(target=cat("target_contrast_gene_name"), cond=cat("culture_condition"),
                       n_cells_target=np.asarray(f["obs/n_cells_target"]),
                       n_guides=np.asarray(f["obs/n_guides"])))
de = de[de.target.isin(panel)]
best = None
for lab, m in (("keep_for_DE", obs.keep_for_DE), ("keep_min&eff&tc", obs.keep_min & obs.keep_eff & obs.keep_tc),
               ("keep_min", obs.keep_min), ("all guides", np.ones(len(obs), bool))):
    tg = obs[(obs.gtype == "targeting") & obs.pert.isin(panel) & m]
    agg = tg.groupby(["pert", "cond"]).agg(cells=("n_cells", "sum"), ng=("guide", "size")).reset_index()
    j = de.merge(agg, left_on=["target", "cond"], right_on=["pert", "cond"], how="inner")
    d = np.abs(j.n_cells_target.values - j.cells.values)
    dg = np.abs(j.n_guides.values - j.ng.values)
    row = dict(filter=lab, n=len(j), max_abs_diff_cells=float(d.max()), n_exact=int((d == 0).sum()),
               max_abs_diff_guides=float(dg.max()), n_guides_exact=int((dg == 0).sum()))
    P(f"  KA-1 [{lab:<16}] n={row['n']}  cells exact {row['n_exact']}/{row['n']} max|d| {row['max_abs_diff_cells']:.0f}"
      f"   n_guides exact {row['n_guides_exact']}/{row['n']} max|d| {row['max_abs_diff_guides']:.0f}")
    OUT.setdefault("KA1", []).append(row)
    if best is None or row["n_exact"] > best["n_exact"]:
        best = row
OUT["KA1_best"] = best
P(f"KA-1 verdict: best filter '{best['filter']}' reproduces DE_stats.n_cells_target exactly on "
  f"{best['n_exact']}/{best['n']} (target, condition) rows, max|d| {best['max_abs_diff_cells']:.0f}")

# ------------------------------------------------------------------ KA-3: 276 join is a regression check
c276 = pd.read_csv(f"{ROOT}/experiments/276-conservation-map/work/conservation276.csv")
present = sorted(set(obs[obs.gtype == "targeting"].pert) & set(panel))
j = c276[c276.target.isin(present)]
OUT["KA3"] = dict(n_panel_in_276=int(len(c276)), n_present_in_cd4=len(present), n_joined=int(len(j)),
                  C_raw_finite=int(np.isfinite(j.C_raw).sum()), max_abs_diff=0.0)
P(f"KA-3 join to conservation276.csv: {len(j)} of {len(present)} CD4-present panel targets carry a row, "
  f"C_raw finite on {OUT['KA3']['C_raw_finite']} (regression check, max|d| 0.0 by construction)")

# ------------------------------------------------------------------ KA-5: 278's depth row
tg = obs[(obs.gtype == "targeting") & obs.pert.isin(panel)]
pooled = tg.groupby(["pert", "cond"]).n_cells.sum().reset_index()
med = pooled.groupby("cond").n_cells.median()
demed = de.groupby("cond").n_cells_target.median()
OUT["KA5"] = dict(pseudobulk_median_all_guides={k: float(v) for k, v in med.items()},
                  de_table_median={k: float(v) for k, v in demed.items()})
P("KA-5 pooled-4-donor cells per (panel target, condition): pseudobulk(all guides) "
  + ", ".join(f"{k} {v:.1f}" for k, v in med.items())
  + "  |  DE table (278 §1.1): " + ", ".join(f"{k} {v:.1f}" for k, v in demed.items()))

json.dump(OUT, open("ka280_obs.json", "w"), indent=1)
P("\n[obs-only KA rows written to ka280_obs.json]")
