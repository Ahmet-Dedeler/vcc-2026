"""813: target coverage of the IGVF TeloHAEC genome-wide CRISPRi release (resting rep1, IGVFFI2522NCFR) from its
guide assignment only (uns/guide_names + obsm/guide_assignment/indices, HTTP range reads, byte-capped; no expression read).

Usage (from this dir): uv run python -I cov813.py <out.json>
Counts are cell-guide incidences: a cell carrying k guides counts once for each of its k targets (MOI ~5 design).
"""
import collections
import csv
import json
import sys

import h5py
import numpy as np

R = "~/Code/virtual-cell-challenge-2026/"
URL = "https://api.data.igvf.org/matrix-files/IGVFFI2522NCFR/@@download/IGVFFI2522NCFR.h5ad"
OUT = sys.argv[1]
src = open(R + "experiments/813-public-data/work/h5peek813.py").read().split("rf = RangeFile(URL)")[0]
ns = {"__name__": "x"}
sys.argv = ["x", "u", "o", "70"]
exec(src, ns)
rf = ns["RangeFile"](URL)
with h5py.File(rf, "r") as f:
    names = [x.decode() if isinstance(x, bytes) else str(x) for x in f["uns/guide_names"][:]]
    idx = f["obsm/guide_assignment/indices"][:]


def tgt(n):
    if n.startswith("non-targeting"):
        return "NTC"
    p = n.split("_")
    return p[1] if len(p) >= 3 else n


T = [tgt(n) for n in names]
per_guide = np.bincount(idx, minlength=len(names))
cnt = collections.Counter()
for t, c in zip(T, per_guide):
    cnt[t] += int(c)
g_per = collections.Counter(T)
genes = [t for t in cnt if t != "NTC"]
v = np.array([cnt[t] for t in genes])
panel = set(r["target_gene"] for r in csv.DictReader(open(R + "data/challenge/pert_counts.csv")))
ess = set(open(R + "data/donors/k562_essential/perts.txt").read().split())
G = set(genes)
out = {"source": URL, "n_guides": len(names), "n_ntc_guides": int(g_per["NTC"]), "n_target_genes": len(genes),
       "guides_per_gene_median": float(np.median([g_per[t] for t in genes])),
       "incidences_total": int(idx.size),
       "cells_per_gene_rep1_quantiles_0_10_25_50_75_90": np.quantile(v, [0, .1, .25, .5, .75, .9]).tolist(),
       "ntc_incidences_rep1": cnt["NTC"],
       "panel_in_library": len(panel & G), "panel_cells_rep1_median": float(np.median([cnt[t] for t in panel & G])),
       "panel_ge100_rep1": int(sum(cnt[t] >= 100 for t in panel & G)), "essential_in_library": len(ess & G),
       "bytes_fetched": rf.fetched}
for k in range(1, 7):
    P = set(r["target_gene"] for r in csv.DictReader(open(R + f"experiments/801-lab-folds/work/L{k}/bundle/pert_counts.csv")))
    out[f"L{k}_in_library"] = len(P & G)
    out[f"L{k}_ge100_rep1"] = int(sum(cnt[t] >= 100 for t in P & G))
json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps(out, indent=1))
