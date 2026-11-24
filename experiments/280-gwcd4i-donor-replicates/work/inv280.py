"""Inventory from obs alone (no X read): guides and cells per (panel target, donor, condition)."""
import numpy as np, pandas as pd, json
ROOT = "~/Code/virtual-cell-challenge-2026"
z = np.load("obs280.npz", allow_pickle=True)
df = pd.DataFrame({k: z[k] for k in ("donor","cond","gtype","run","pert","guide","n_cells","total",
                                     "keep_for_DE","keep_min","keep_eff","keep_tc","keep_tg")})
panel = sorted(set(l.split(",")[0] for l in open(f"{ROOT}/data/challenge/pert_counts.csv").read().strip().split("\n")[1:]))
DON = {"CE0008162":"D1","CE0010866":"D2","CE0008678":"D3","CE0006864":"D4"}
df["D"] = df.donor.map(DON)
tg = df[(df.gtype=="targeting") & df.pert.isin(panel)]
print("panel targeting rows:", len(tg), "| keep_for_DE:", int(tg.keep_for_DE.sum()))
print("\n-- keep flag breakdown on panel targeting rows --")
for c in ("keep_min","keep_eff","keep_tc","keep_tg","keep_for_DE"):
    print(f"  {c:<16} {int(tg[c].sum()):>6} / {len(tg)}")
for flt, lab in ((tg, "ALL guides"), (tg[tg.keep_for_DE], "keep_for_DE")):
    g = flt.groupby(["cond","D","pert"]).agg(n_guides=("guide","size"), cells=("n_cells","sum"))
    print(f"\n== {lab}: guides per (target, donor, condition) ==")
    print(g.n_guides.value_counts().sort_index().to_string())
    print(f"   cells per (target,donor,cond): median {g.cells.median():.0f}  q10 {g.cells.quantile(.1):.0f}  q90 {g.cells.quantile(.9):.0f}  >=100: {(g.cells>=100).mean():.3f}")
    t = flt.groupby(["cond","pert"]).agg(cells=("n_cells","sum"), nd=("D","nunique"))
    print(f"   pooled-4-donor cells per (target,cond): median {t.cells.median():.0f} ; targets with all 4 donors: {(t.nd==4).sum()} / {len(t)}")
ntc = df[df.gtype!="targeting"]
n = ntc.groupby(["cond","D"]).agg(guides=("guide","size"), cells=("n_cells","sum"))
print("\n== non-targeting pool per (donor, condition) ==\n", n.to_string())
print("   NTC cells per guide: median %.0f q10 %.0f q90 %.0f" % (ntc.n_cells.median(), ntc.n_cells.quantile(.1), ntc.n_cells.quantile(.9)))
miss = sorted(set(panel) - set(df.pert.unique()))
print("\npanel targets absent from pseudobulk:", miss)
