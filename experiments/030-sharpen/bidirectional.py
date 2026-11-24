"""030f — the third basis: reweight along the PERTURBATION axis too.

Two smooth reweightings have paid, in two different bases:

  directions across perturbations   spectral whitening      +0.070
  magnitudes within a perturbation  power sharpening        +0.038

There is one obvious basis left: **a gene's profile across the panel**. For gene g, the 300
values `e_{p,g}` say how much that gene responds to each perturbation, and the same argument
applies — the large ones are the better-estimated ones. `gene_norm` already touches this axis,
but only as a diagonal rescaling by each gene's mean |effect|; it does not reshape the
distribution.

Doing both axes at once is Sinkhorn-style bi-normalisation, which is standard in retrieval
(and is what TF-IDF approximates in one step). Three variants:

  row     sharpen within each perturbation                        (the incumbent lever)
  col     sharpen each gene's profile across the panel
  both    alternate, Sinkhorn-style, to a fixed point

If the column axis pays as well, the mechanism generalises and there may be more of these. If
it does not, that is a real boundary on where the redundancy lives -- and worth knowing before
looking for a fourth.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments" / "018-proxy6"))
from proxy6 import pds_cosine   # noqa: E402
import run as proxy_run          # noqa: E402

GAMMA = 0.65


def sharpen(M, q):
    return M if q == 1.0 else np.sign(M) * np.abs(M) ** q


def col_sharpen(M, r):
    """Reshape each gene's across-panel profile, holding its column scale fixed."""
    if r == 1.0:
        return M
    A = np.abs(M)
    s0 = A.mean(0) + 1e-12
    out = np.sign(M) * (A ** r)
    return out * (s0 / (np.abs(out).mean(0) + 1e-12))


folds, common = proxy_run.build_folds("1:1")
print(f"panel {len(common)}  gamma={GAMMA}\n")
res = {}
for h, f in folds.items():
    P0 = f["W"].apply(f["P"], GAMMA)
    T = f["T"]
    for q in (1.0, 1.5, 2.0):
        for r in (1.0, 1.25, 1.5, 2.0):
            res.setdefault(f"row q={q:g} col r={r:g}", {})[h] = \
                pds_cosine(col_sharpen(sharpen(P0, q), r), T)
    # Sinkhorn-style: alternate mild sharpening on both axes to a fixed point
    for n_it in (2, 4):
        M = P0.copy()
        for _ in range(n_it):
            M = sharpen(M, 1.35)
            M = col_sharpen(M, 1.35)
        res.setdefault(f"sinkhorn x{n_it} (1.35)", {})[h] = pds_cosine(M, T)

print(f"{'config':24} " + "  ".join(f"{k.replace('orion_','')[:9]:>10}" for k in folds)
      + f" {'macro':>8}")
out = {}
for k, v in sorted(res.items()):
    macro = float(np.mean(list(v.values())))
    out[k] = {**v, "macro": macro}
    print(f"{k:24} " + "  ".join(f"{v[n]:10.4f}" for n in folds) + f" {macro:8.4f}")
best = max(out, key=lambda k: out[k]["macro"])
print(f"\nbest: {best} -> {out[best]['macro']:.4f}")
print(f"row-only reference (q=2, r=1): {out['row q=2 col r=1']['macro']:.4f}")
Path("experiments/030-sharpen/bidirectional.json").write_text(json.dumps(out, indent=2))
