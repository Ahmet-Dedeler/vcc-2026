"""282 metadata probe -- design counts only.  Reads obs280.npz (.obs of the pseudobulk) and the HGNC
table.  NO expression value is opened; nothing here can see a state-dependence number."""
from __future__ import annotations
import csv, json
import numpy as np, pandas as pd

ROOT = "~/Code/virtual-cell-challenge-2026"
W280 = f"{ROOT}/experiments/280-gwcd4i-donor-replicates/work"
DONORS = {"CE0008162": "D1", "CE0010866": "D2", "CE0008678": "D3", "CE0006864": "D4"}
DORDER = ["D1", "D2", "D3", "D4"]
CONDS = ["Rest", "Stim8hr", "Stim48hr"]

z = np.load(f"{W280}/obs280.npz", allow_pickle=True)
sel = np.load(f"{W280}/bank280_rows.npy")
o = pd.DataFrame({k: z[k][sel] for k in ("donor","cond","gtype","run","pert","guide","n_cells","total","keep_for_DE")})
for c in ("donor","cond","gtype","pert","guide"): o[c] = o[c].astype(str)
o["D"] = o.donor.map(DONORS)
panel = [l.split(",")[0] for l in open(f"{ROOT}/data/challenge/pert_counts.csv").read().strip().split("\n")[1:]]
P = set(panel)
tg = o[(o.gtype == "targeting") & o.pert.isin(P)]
out = {}
for keep in (True, False):
    t2 = tg[tg.keep_for_DE] if keep else tg
    have = {}
    for (c, d, t), g in t2.groupby(["cond","D","pert"]):
        have.setdefault(t, set()).add((c, d))
    need = {(c, d) for c in CONDS for d in DORDER}
    common = sorted(t for t, s in have.items() if need <= s)
    per_cond = {c: sorted(t for t, s in have.items() if all((c, d) in s for d in DORDER)) for c in CONDS}
    out["keepDE" if keep else "allguides"] = dict(
        n_common_all3=len(common), n_per_cond={c: len(v) for c, v in per_cond.items()},
        common=common)
    print(("keep_for_DE" if keep else "all guides"), "targets in all 4 donors x all 3 states:", len(common),
          " per-cond (280's sets):", {c: len(v) for c, v in per_cond.items()}, flush=True)

common = out["keepDE"]["common"]
# cells per (target, state) pooled over 4 donors, and per donor-PAIR
t2 = tg[tg.keep_for_DE & tg.pert.isin(common)]
for c in CONDS:
    g = t2[t2.cond == c].groupby("pert").n_cells.sum()
    gp = t2[(t2.cond == c) & t2.D.isin(["D1","D2"])].groupby("pert").n_cells.sum()
    print(f"  {c:<9} cells/(target,4 donors) median {g.median():.0f} q10 {g.quantile(.1):.0f} q90 {g.quantile(.9):.0f}"
          f" | per D1+D2 pair median {gp.median():.0f}", flush=True)
# guide-set stability across states
same = 0
for t in common:
    sets = {c: {(d, gg) for d, gg in t2[(t2.cond == c) & (t2.pert == t)][["D","guide"]].values} for c in CONDS}
    if sets["Rest"] == sets["Stim8hr"] == sets["Stim48hr"]:
        same += 1
print(f"  targets whose (donor,guide) keep_for_DE set is IDENTICAL in all 3 states: {same}/{len(common)}", flush=True)

# HGNC families of the panel
fam = {}
with open(f"{ROOT}/data/external/xcell/hgnc_protein_coding.txt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        fam[r["symbol"]] = ([x.strip() for x in (r["gene_group"] or "").split("|") if x.strip()], r["name"], r["location"])
n_fam = {}
for t in common:
    for f in fam.get(t, ([], "", ""))[0]:
        n_fam[f] = n_fam.get(f, 0) + 1
big = sorted(n_fam.items(), key=lambda kv: -kv[1])[:22]
print("  HGNC families among the common targets (>=5):", [(k, v) for k, v in big if v >= 5], flush=True)
print("  common targets with NO hgnc row:", sum(1 for t in common if t not in fam), flush=True)
json.dump({k: (v if k != "keepDE" and k != "allguides" else {kk: vv for kk, vv in v.items()})
           for k, v in out.items()}, open("meta282.json", "w"), indent=1)
