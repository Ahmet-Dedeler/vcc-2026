"""024b — the covariance was under-sampled, so the gamma optimum was under-estimated.

Capping the basis at 6,000 perturbations per donor was arbitrary. The donors carry 9,787 /
17,767 / 17,977. Going to the full catalogue moves macro PDS 0.7995 -> 0.8071 and pushes the
optimal gamma from 0.65 to at least 0.75 -- which is the expected coupling: a better-
conditioned covariance can be inverted harder before noise in the small eigenvalues starts to
dominate. The 0.75 point was the edge of the previous grid, so the optimum needs re-finding.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.donor_transfer import Donor, _norm_scale   # noqa: E402
from vcc2026.whiten import Whitener                      # noqa: E402

DONORS = {"k562_gwps": "data/donors/k562_gwps",
          "orion_hct116": "data/donors/orion_hct116",
          "orion_hek293t": "data/donors/orion_hek293t"}
GAMMAS = [0.65, 0.75, 0.85, 0.95, 1.05, 1.2]


def pds(pred, true):
    P = pred / (np.linalg.norm(pred, axis=1, keepdims=True) + 1e-12)
    T = true / (np.linalg.norm(true, axis=1, keepdims=True) + 1e-12)
    D = 1.0 - P @ T.T
    n = D.shape[1]
    sd = np.diag(D)[:, None]
    return float((1.0 - ((D < sd).sum(1) + 0.5 * ((D == sd).sum(1) - 1)) / (n - 1)).mean())


banks = {}
for name, path in DONORS.items():
    d = Donor.load(path)
    _, s = d.decompose()
    gl = [str(g) for g in d.genes]
    banks[name] = dict(s=(s / _norm_scale(s)).astype(np.float32), raw=s,
                       gpos={g: i for i, g in enumerate(gl)}, gset=set(gl),
                       idx={p: i for i, p in enumerate(d.perts)})
vcc = list(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])
common = [g for g in vcc if all(g in b["idx"] for b in banks.values())]
print(f"panel {len(common)}\n")

res = {}
for held in banks:
    others = [n for n in banks if n != held]
    shared = sorted(set.intersection(*[banks[n]["gset"] for n in [held] + others]))
    feat = [g for g in shared if g not in set(common)]
    col = {n: np.fromiter((banks[n]["gpos"][g] for g in feat), np.int64, len(feat))
           for n in [held] + others}
    T = np.vstack([banks[held]["raw"][banks[held]["idx"][p]][col[held]]
                   for p in common]).astype(np.float64)
    sc = {n: float(np.median(np.linalg.norm(banks[n]["s"][:, col[n]], axis=1))) or 1.0
          for n in others}
    ref = float(np.median(list(sc.values())))
    P = sum((ref / sc[n]) * np.vstack(
        [banks[n]["s"][banks[n]["idx"][p]][col[n]] for p in common]).astype(np.float64)
        for n in others) / len(others)
    B = np.vstack([(ref / sc[n]) * banks[n]["s"][:, col[n]].astype(np.float64)
                   for n in others])
    print(f"fold {held:14} full basis {B.shape}", flush=True)
    W = Whitener.fit(B)
    del B
    for g in GAMMAS:
        res.setdefault(f"full g={g:g}", {})[held] = pds(W.apply(P, g), T)

print(f"\n{'config':16} " + "  ".join(f"{k.replace('orion_','')[:9]:>10}" for k in banks)
      + f" {'macro':>8}")
out = {}
for k, v in sorted(res.items()):
    macro = float(np.mean(list(v.values())))
    out[k] = {**v, "macro": macro}
    print(f"{k:16} " + "  ".join(f"{v[n]:10.4f}" for n in banks) + f" {macro:8.4f}")
best = max(out, key=lambda k: out[k]["macro"])
print(f"\nbest: {best} -> {out[best]['macro']:.4f}")
Path("experiments/024-whiten-variants/bigbasis.json").write_text(json.dumps(out, indent=2))
