"""036 — what is the marginal value of one more donor line, under whitening?

Every closed avenue now points the same way: the binding constraint is directional accuracy,
and the only two things that have ever moved it are **more donor lines** and **better
geometry**. The geometry family is exhausted (exp 030f). So the question is whether acquiring a
fourth line is worth the day it would cost.

There is one candidate within reach — Feng et al. 2025, a genome-scale CRISPRi atlas across 34
iPSC lines (7,226 genes targeted, ~2M cells, Figshare). It is *pluripotent* where the contexts
are somatic, which is the same mismatch that made H1 a poor donor, and Bio_NUDT runs
`gwps+feng` at `pds` 0.665. So the prior is not good and the cost is high.

Rather than guess, measure the marginal value of the Nth donor directly. The 4-line LODO gives
folds with **three** predictors available, so 1 -> 2 -> 3 can be measured under whitening on the
same folds, and the curve extrapolated to 4.

Experiment 019 already showed donors are worth *more* once whitening is on (the second donor
went from +0.017 to +0.033), because averaging's cost was the smoothing whitening removes. If
the curve is still steep at 3, Feng is worth chasing despite the tissue mismatch. If it has
flattened, it is not.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.donor_transfer import Donor, _norm_scale   # noqa: E402
from vcc2026.whiten import Whitener                      # noqa: E402

LINES = {"k562_gwps": "data/donors/k562_gwps", "orion_hct116": "data/donors/orion_hct116",
         "orion_hek293t": "data/donors/orion_hek293t", "h1_2025": "data/donors/h1_2025"}
GAMMA = 0.65

banks = {}
for n, p in LINES.items():
    d = Donor.load(p)
    _, s = d.decompose()
    gl = [str(g) for g in d.genes]
    banks[n] = dict(s=(s / _norm_scale(s)).astype(np.float32), raw=s,
                    gpos={g: i for i, g in enumerate(gl)}, gset=set(gl),
                    idx={q: i for i, q in enumerate(d.perts)})
common = sorted(set.intersection(*[set(b["idx"]) for b in banks.values()]))
print(f"{len(common)} perturbations across all 4 lines\n")


def pds(P, T):
    Pn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12)
    Tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    D = 1.0 - Pn @ Tn.T
    n = D.shape[1]
    sd = np.diag(D)[:, None]
    return float((1.0 - ((D < sd).sum(1) + 0.5 * ((D == sd).sum(1) - 1)) / (n - 1)).mean())


res = {1: [], 2: [], 3: []}
for held in banks:
    others = [n for n in banks if n != held]
    shared = sorted(set.intersection(*[banks[n]["gset"] for n in [held] + others]))
    feat = [g for g in shared if g not in set(common)]
    col = {n: np.fromiter((banks[n]["gpos"][g] for g in feat), np.int64, len(feat))
           for n in [held] + others}
    T = np.vstack([banks[held]["raw"][banks[held]["idx"][q]][col[held]]
                   for q in common]).astype(np.float64)
    sc = {n: float(np.median(np.linalg.norm(banks[n]["s"][:, col[n]], axis=1))) or 1.0
          for n in others}
    ref = float(np.median(list(sc.values())))
    parts = {n: (ref / sc[n]) * np.vstack(
        [banks[n]["s"][banks[n]["idx"][q]][col[n]] for q in common]).astype(np.float64)
        for n in others}
    full = {n: (ref / sc[n]) * banks[n]["s"][:, col[n]].astype(np.float64) for n in others}
    print(f"fold {held:14} fitting whiteners", flush=True)
    # One whitener per predictor SUBSET: the operator is estimated from exactly the lines the
    # prediction uses, so a k-donor configuration is never given a k+1-donor covariance.
    for k in (1, 2, 3):
        for combo in itertools.combinations(others, k):
            W = Whitener.fit(np.vstack([full[n] for n in combo]))
            P = sum(parts[n] for n in combo) / k
            res[k].append(pds(W.apply(P, GAMMA), T))

print(f"\n{'donors':>7} {'folds x combos':>15} {'macro PDS':>11} {'marginal':>10}")
prev = None
out = {}
for k in (1, 2, 3):
    m = float(np.mean(res[k]))
    out[k] = dict(mean=m, n=len(res[k]), sd=float(np.std(res[k])))
    print(f"{k:7} {len(res[k]):15} {m:11.4f} "
          f"{'' if prev is None else f'{m - prev:+10.4f}'}")
    prev = m
d21 = out[2]["mean"] - out[1]["mean"]
d32 = out[3]["mean"] - out[2]["mean"]
ratio = d32 / d21 if d21 else float("nan")
print(f"\ndecay ratio (3rd donor / 2nd donor) = {ratio:.2f}")
print(f"extrapolated 4th donor              = {d32 * ratio:+.4f} macro PDS")
print(f"for scale: whitening itself is worth +0.070, and the 2nd donor here is {d21:+.4f}")
Path("experiments/036-donor-marginal/results.json").write_text(json.dumps(out, indent=2))
