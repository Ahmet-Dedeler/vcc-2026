"""436 (c) side question: do the essential screens carry reproducible direction on the targets an ESSENTIAL-FREE panel
could contain (targets in the screens' union that are NOT on the K562-essential list)? Those are the only targets they
can cover in October if D/E/F follows Arc's criterion at all loosely (238's efp/efr: 9 each, all off-list).
Cross-screen reproducibility: for a target in >= 2 screens, for each ordered pair (A, B) the sign agreement of B with A
on A's top-200 |delta| genes that both measure (0.5 = no shared direction). Split on-list / off-list; also the
effect-row norm. Writes work/xscreen436.json."""
from __future__ import annotations

import csv
import json
from itertools import permutations
from pathlib import Path

import numpy as np

import signs436 as S6

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ess_list = set(S6._lines(S6.ESS_DIRS["k562_essential"] / "perts.txt"))
tabs = {}
for nm, d in S6.ESS_DIRS.items():
    P = S6._lines(d / "perts.txt"); G = S6._lines(d / "genes.txt")
    tabs[nm] = ({p: i for i, p in enumerate(P)}, G, np.load(d / "delta.npy", mmap_mode="r"))
common = sorted(set.intersection(*[set(t[1]) for t in tabs.values()]))
cidx = {nm: np.array([{g: i for i, g in enumerate(t[1])}[g] for g in common]) for nm, t in tabs.items()}
union = sorted(set().union(*[set(t[0]) for t in tabs.values()]))


def read_col(p, col):
    with open(p) as f:
        return [row[col] for row in csv.DictReader(f)]


efcov = set()
for b in ("ef-bundle", "ef-bundle-r"):
    efcov |= set(read_col(ROOT / f"experiments/238-essential-free-panel/work/{b}/pert_counts.csv", "target_gene")) & set(union)
res = {"on": [], "off": [], "ef": []}
norms = {"on": [], "off": []}
for t in union:
    rows = {nm: np.asarray(tb[2][tb[0][t]], dtype=np.float64)[cidx[nm]] for nm, tb in tabs.items() if t in tb[0]}
    grp = "on" if t in ess_list else "off"
    for r_ in rows.values():
        norms[grp].append(float(np.linalg.norm(r_)))
    if len(rows) < 2:
        continue
    ag = []
    for a, b in permutations(rows, 2):
        top = np.argsort(-np.abs(rows[a]))[:200]
        ag.append(float(np.mean(np.sign(rows[a][top]) == np.sign(rows[b][top]))))
    res[grp].append(float(np.mean(ag)))
    if t in efcov:
        res["ef"].append((t, float(np.mean(ag)), len(rows)))
out = {"common_genes": len(common), "union": len(union)}
for g in ("on", "off"):
    v = np.array(res[g])
    out[g] = {"n_targets_ge2_screens": int(v.size), "agree_median": float(np.median(v)), "agree_q25": float(np.quantile(v, .25)),
              "agree_q75": float(np.quantile(v, .75)), "frac_agree_ge_0.7": float(np.mean(v >= 0.7)),
              "row_norm_median": float(np.median(norms[g]))}
    print(f"{g:3s}-list targets in >= 2 screens: {v.size}; top-200 cross-screen sign agreement median {np.median(v):.3f} "
          f"(IQR {np.quantile(v, .25):.3f}-{np.quantile(v, .75):.3f}), >= 0.7 in {np.mean(v >= 0.7):.2f}; row norm median {np.median(norms[g]):.2f}")
out["ef_covered"] = res["ef"]
print("238 efp/efr covered targets (all off-list):", [(t, round(a, 3), n) for t, a, n in res["ef"]])
(HERE / "work/xscreen436.json").write_text(json.dumps(out, indent=1))
