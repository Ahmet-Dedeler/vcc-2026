"""480 resume-2: the truth's expression-direction profile per surface (the Orion "technical blob" check of §5.1).

For each surface's scorer-rule gate: the share of gate cells (target x gene) whose truth lfc is up, by the gene's
control expression (count per cell = control CPM x the surface's median library / 1e6), in deciles and tertiles of the
gate-pooled distribution. H1 leans down with expression (0.51 -> 0.41, extras480 updown); a surface whose low-expressed
gate genes are mostly up and high-expressed mostly down carries a target-independent artefact. The thinned H1 gates
(h1thin480.py levels: cells per target x depth fraction) separate a detection-power asymmetry (decreases of
low-expressed genes are hard to call with few cells / little depth) from a line- or study-specific blob.

    uv run python artefact480.py     # -> work/artefact480.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def surface(name):
    if name == "h1":
        r = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
        keep = r["keep"].astype(bool)
        genes, targets = r["genes"].astype(str)[keep], r["targets"].astype(str)
        q, lfc = r["q"].astype(float), r["lfc"].astype(float)
        ref_lin = r["ref_lin"].astype(float)
        cnt = ref_lin * float(np.median(r["libraries"].astype(float))) / 1e6
    else:
        z = np.load(HERE / f"work/gate-{name}.npz", allow_pickle=True)
        genes, targets = z["genes"].astype(str), z["targets"].astype(str)
        q, lfc = z["q"].astype(float), z["lfc"].astype(float)
        cnt = z["ctrl_cpm"].astype(float) * float(z["lib_med"]) / 1e6
    gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & (genes[None, :] != targets[:, None])
    gate &= (gate.sum(1) >= 10)[:, None]
    return gate, lfc, np.broadcast_to(cnt, lfc.shape)


out = {}
for name in ["h1", "h1_full_ka", "h1_c250", "h1_c400_d56", "h1_c150", "h1_c250_d56", "h1_c100", "h1_thin",
             "h1_c150_d56", "orion_hct116", "orion_hek293t", "k562_gwps"]:
    gate, lfc, C = surface(name)
    up, c = (lfc[gate] > 0), C[gate]
    res = {}
    for nb, key in [(10, "decile"), (3, "tertile")]:
        edges = np.quantile(c, np.linspace(0, 1, nb + 1)[1:-1])
        b = np.digitize(c, edges)
        res[key] = [float(up[b == i].mean()) for i in range(nb)]
    res["all"] = float(up.mean())
    out[name] = res
    print(f"{name:14s} P(up) all {res['all']:.3f} | tertiles " + " ".join(f"{v:.3f}" for v in res["tertile"]) +
          " | deciles " + " ".join(f"{v:.2f}" for v in res["decile"]))
json.dump(out, open(HERE / "work/artefact480.json", "w"), indent=1)
