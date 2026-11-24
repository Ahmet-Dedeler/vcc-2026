"""Summarise de_power_panel_targets.csv: by cell bin, on the 50 candidates, gate>=10 counts, curve."""
import json, numpy as np, pandas as pd
from pathlib import Path
W = Path("~/Code/virtual-cell-challenge-2026/experiments/260-raw-de-diagnostic-feasibility/work")
df = pd.read_csv(W / "de_power_panel_targets.csv")
cov = pd.read_csv(W / "panel_targets_coverage.csv")
cand = pd.read_csv(W / "candidate_targets_ge200.csv")
df = df.merge(cov[["target", "six_union", "orion_both", "cd4_complete", "h1_pair", "pm"]], on="target", how="left")
print(f"targets scored: {len(df)}; scorer convention = tie-corrected (cell_eval2/de_compute.py:871-872 scanpy wilcoxon tie_correct=True)")
print("control median UMI", json.load(open(W / "de_power_summary.json"))["ctrl_median_umi"])
rows = []
for lo, hi in ((0, 100), (100, 150), (150, 200), (200, 250), (250, 300), (300, 401)):
    d = df[(df.cells >= lo) & (df.cells < hi)]
    if len(d) == 0: continue
    rows.append(dict(bin=f"[{lo},{hi})", n=len(d), n_real_median=d.n_real_tc.median(), q1=d.n_real_tc.quantile(.25), q3=d.n_real_tc.quantile(.75),
                     gate10=int((d.gate_tc >= 10).sum()), gate10_frac=(d.gate_tc >= 10).mean(), n_real_uc_median=d.n_real_uc.median(),
                     frac_neg=d.frac_neg_tc.median()))
t = pd.DataFrame(rows); print(t.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
print("\nALL 272: gate>=10 (tc):", int((df.gate_tc >= 10).sum()), "median n_real", df.n_real_tc.median(), "frac_neg median", round(df.frac_neg_tc.median(), 3),
      "zero-DE targets:", int((df.n_real_tc == 0).sum()))
c = df[df.target.isin(cand.target)]
print(f"\n50 CANDIDATES (>=200 cells, union+Orion+CD4): scored {len(c)}; n_real tc median {c.n_real_tc.median():.0f} IQR [{c.n_real_tc.quantile(.25):.0f},{c.n_real_tc.quantile(.75):.0f}] min {c.n_real_tc.min()} max {c.n_real_tc.max()}; gate>=10: {(c.gate_tc>=10).sum()}/{len(c)}; >=30: {(c.n_real_tc>=30).sum()}; >=89: {(c.n_real_tc>=89).sum()}")
print("  candidates below the gate:", ", ".join(f"{r.target}({r.cells}c,{r.n_real_tc})" for r in c[c.gate_tc < 10].itertuples()))
print("  candidates sorted by n_real:", ", ".join(f"{r.target}:{r.n_real_tc}" for r in c.sort_values('n_real_tc', ascending=False).itertuples()))
c2 = df[(df.cells >= 150) & df.six_union & df.orion_both & df.cd4_complete]
print(f"\n91 at >=150 cells (union+Orion+CD4): gate>=10: {(c2.gate_tc>=10).sum()}/{len(c2)}; median n_real {c2.n_real_tc.median():.0f}")
k = df[df.target_in_axis & df.target_knockdown_lfc.notna()]
print(f"\ntarget-gene knockdown (targets on the axis): n={len(k)}, median lfc {k.target_knockdown_lfc.median():.2f}, frac < -0.5: {(k.target_knockdown_lfc < -0.5).mean():.2f}, frac < 0: {(k.target_knockdown_lfc < 0).mean():.2f}")
cur = pd.read_csv(W / "de_power_curve.csv"); print("\nsubsampling curve (n_real tc vs cells):"); print(cur.pivot(index="target", columns="cells", values="n_real_tc").to_string())
print("\nchallenge reference points: h1-test median n_real 89 at 400 cells / 20k UMI (188); challenge nmae cohort 209-261/300 (metrics brief §8)")
c.sort_values("n_real_tc", ascending=False).to_csv(W / "candidates_de_power.csv", index=False)
