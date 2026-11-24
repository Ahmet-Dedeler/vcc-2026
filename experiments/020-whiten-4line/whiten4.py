"""020 — does whitening still pay with three predictors?

Experiments 017/019 were run on the 3-line LODO, which leaves **two** predictors per fold.
The submission uses **three** donors. Two things could differ:

1. Averaging three sources cancels more independent noise than averaging two, so the
   prediction fed to the whitener is already less noisy and the operator has less to fix.
2. Whitening's whole mechanism is removing the shared component that averaging amplifies.
   With three sources that shared component is *more* dominant, so it could matter more.

The 4-line LODO built for experiment 016 (K562-gwps, Orion-HCT116, Orion-HEK293T, H1-2025)
gives three predictors per fold. Its caveat carries over unchanged and is important: of the
226 perturbations shared by all four lines, only 20 are on the VCC panel, because H1-2025's
panel is essential-gene enriched and the VCC panel excludes exactly that class. So this
measures whitening on a *different gene population* than the task. It is a robustness check
on the mechanism, not a second estimate of the effect size.
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

LINES = {"k562_gwps": "data/donors/k562_gwps", "orion_hct116": "data/donors/orion_hct116",
         "orion_hek293t": "data/donors/orion_hek293t", "h1_2025": "data/donors/h1_2025"}

banks = {}
for n, p in LINES.items():
    d = Donor.load(p)
    _, s = d.decompose()
    gl = [str(g) for g in d.genes]
    banks[n] = dict(s=(s / _norm_scale(s)).astype(np.float32), raw=s,
                    gpos={g: i for i, g in enumerate(gl)}, gset=set(gl),
                    idx={q: i for i, q in enumerate(d.perts)})

common = sorted(set.intersection(*[set(b["idx"]) for b in banks.values()]))
vcc = set(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])
print(f"{len(common)} perturbations across all 4 lines ({len(set(common) & vcc)} on the "
      f"VCC panel)\n")


def pds(pred, true):
    P = pred / (np.linalg.norm(pred, axis=1, keepdims=True) + 1e-12)
    T = true / (np.linalg.norm(true, axis=1, keepdims=True) + 1e-12)
    D = 1.0 - P @ T.T
    n = D.shape[1]
    sd = np.diag(D)[:, None]
    return float((1.0 - ((D < sd).sum(1) + 0.5 * ((D == sd).sum(1) - 1)) / (n - 1)).mean())


GAMMAS = [0.0, 0.4, 0.55, 0.65, 0.8]
WEIGHTS = {"equal": [1, 1, 1], "primary 4:1:1": [4, 1, 1]}
rows = {}
rng = np.random.default_rng(0)
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
    parts, Bl = {}, []
    for n in others:
        parts[n] = (ref / sc[n]) * np.vstack(
            [banks[n]["s"][banks[n]["idx"][q]][col[n]] for q in common]).astype(np.float64)
        M = banks[n]["s"][:, col[n]]
        if M.shape[0] > 6000:
            M = M[np.sort(rng.choice(M.shape[0], 6000, False))]
        Bl.append((ref / sc[n]) * M.astype(np.float64))
    print(f"fold {held:14} {T.shape[0]} perts x {T.shape[1]} genes, "
          f"whitener basis {sum(b.shape[0] for b in Bl)}", flush=True)
    W = Whitener.fit(np.vstack(Bl))
    for wname, ws in WEIGHTS.items():
        w = dict(zip(others, ws))
        P = sum((w[n] / sum(w.values())) * parts[n] for n in others)
        for g in GAMMAS:
            rows.setdefault(f"{wname} gamma={g:g}", {})[held] = pds(W.apply(P, g), T)

print(f"\n{'config':24} " + "  ".join(f"{k.replace('orion_','')[:9]:>10}" for k in banks)
      + f" {'macro':>8}")
out = {}
for k, v in rows.items():
    macro = float(np.mean(list(v.values())))
    out[k] = {**v, "macro": macro}
    print(f"{k:24} " + "  ".join(f"{v[n]:10.4f}" for n in banks) + f" {macro:8.4f}")
best = max(out, key=lambda k: out[k]["macro"])
print(f"\nbest: {best} -> {out[best]['macro']:.4f}")
Path("experiments/020-whiten-4line/results.json").write_text(json.dumps(out, indent=2))
