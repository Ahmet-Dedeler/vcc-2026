"""Does whitening raise the transfer COSINE, or only the ranking?

`pds` is a rank; `expr_mse` is not. The whole "mse is structurally unreachable" argument rests
on one number: `min_alpha mse = 1 - cos^2`, which beats the 0.989 baseline only above
cos = 0.105, and the measured donor-transfer cosine was ~0.08. If whitening improved the
ranking by reshaping the space *without* improving alignment, that argument stands unchanged.
If it also raised the cosine, a sixth metric comes back into play.

Measured in the space the scorer actually reads -- raw log-pseudobulk effects on the held-out
line -- so it is the quantity the mse formula wants, not a whitened-space artefact.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, _norm_scale   # noqa: E402
from vcc2026.whiten import Whitener                      # noqa: E402

DONORS = {"k562_gwps": "data/donors/k562_gwps",
          "orion_hct116": "data/donors/orion_hct116",
          "orion_hek293t": "data/donors/orion_hek293t"}

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


def rowcos(P, T):
    Pn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12)
    Tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    return (Pn * Tn).sum(1)


print(f"{'fold':16} {'cos g=0':>9} {'cos g=0.85':>11} {'min mse g=0':>12} "
      f"{'min mse g=.85':>14}")
rows = []
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
    W = Whitener.fit(np.vstack([(ref / sc[n]) * banks[n]["s"][:, col[n]].astype(np.float64)
                                for n in others]))
    c0, c1 = rowcos(P, T).mean(), rowcos(W.apply(P, 0.85), T).mean()
    rows.append((c0, c1))
    print(f"{held:16} {c0:9.4f} {c1:11.4f} {1-c0**2:12.4f} {1-c1**2:14.4f}", flush=True)

c0 = float(np.mean([r[0] for r in rows])); c1 = float(np.mean([r[1] for r in rows]))
print(f"\nmacro cosine {c0:.4f} -> {c1:.4f}   ({(c1/c0-1)*100:+.1f}%)")
print(f"min mse      {1-c0**2:.4f} -> {1-c1**2:.4f}   (official baseline b = 0.989; "
      f"scaled mse is positive only below it)")
print(f"cos needed for a positive scaled mse: 0.105")
