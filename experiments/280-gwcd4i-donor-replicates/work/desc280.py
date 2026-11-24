"""280 §5 -- DESCRIPTIVE, POST-HOC, no pass rule: which panel targets are donor-variable, and whether donor
variability tracks 276's low-conservation tail.  Also the cross-check against the authors' own
`donor_correlation_*` columns in the on-disk DE table."""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = "~/Code/virtual-cell-challenge-2026"
TAG = "keepDE"
P = lambda *a: print(*a, flush=True)  # noqa: E731

c276 = pd.read_csv(f"{ROOT}/experiments/276-conservation-map/work/conservation276.csv").set_index("target")
de = pd.read_csv("de_stats_panel280.csv")
rows = []
for cond in ("Rest", "Stim8hr", "Stim48hr"):
    z = np.load(f"q1_pertarget_{cond}_{TAG}.npz", allow_pickle=True)
    tt = [str(x) for x in z["targets"]]
    df = pd.DataFrame(dict(target=tt, cond=cond, cos=z["cos"]))
    for k in z.files:
        if "|" in k:
            df[k] = z[k]
    rows.append(df)
D = pd.concat(rows, ignore_index=True)
D = D.merge(c276[["C_raw", "C_rel", "essential", "shared_programme_share", "ccle_median_tpm"]],
            left_on="target", right_index=True, how="left")
D = D.merge(de[["target", "cond", "n_cells_target", "n_total_de_genes", "ontarget_effect_size",
                "donor_correlation_all_mean", "guide_correlation_all"]],
            on=["target", "cond"], how="left")
D.to_csv("desc280_per_target.csv", index=False)

for cond in ("Rest", "Stim8hr", "Stim48hr"):
    d = D[D.cond == cond].dropna(subset=["cos"])
    P(f"\n== {cond}: n={len(d)} ==")
    lo = d.nsmallest(12, "cos")[["target", "cos", "C_raw", "n_cells_target", "n_total_de_genes"]]
    hi = d.nlargest(12, "cos")[["target", "cos", "C_raw", "n_cells_target", "n_total_de_genes"]]
    P("  most donor-VARIABLE (lowest across-donor cosine):")
    P("    " + ", ".join(f"{r.target} {r.cos:+.3f}" for r in lo.itertuples()))
    P("  most donor-REPRODUCIBLE:")
    P("    " + ", ".join(f"{r.target} {r.cos:+.3f}" for r in hi.itertuples()))
    for col in ("C_raw", "C_rel", "n_cells_target", "n_total_de_genes", "shared_programme_share",
                "donor_correlation_all_mean", "guide_correlation_all", "ontarget_effect_size"):
        v = d[col].values.astype(float)
        m = np.isfinite(v)
        if m.sum() < 10:
            continue
        rho, p = spearmanr(v[m], d.cos.values[m])
        P(f"    Spearman(cos, {col:<26}) {rho:+.3f}  p={p:.2g}  n={m.sum()}")
    # does donor variability sit in 276's low-conservation tail?
    v = d.C_raw.values.astype(float)
    m = np.isfinite(v)
    bot30 = v <= np.nanquantile(v, 0.30)
    P(f"    median cos in 276's bottom-30 % C_raw: {np.nanmedian(d.cos.values[bot30 & m]):+.4f}"
      f"   vs the rest {np.nanmedian(d.cos.values[~bot30 & m]):+.4f}   (n {int((bot30&m).sum())}/{int((~bot30&m).sum())})")
