"""026 — re-express donor effects at the target's expression level.

## The mis-specification

Every effect in this repo is `delta = log1p(TS*p_pert) - log1p(TS*p_ctrl)` with TS = 5e4,
which is what the scorer computes. `log1p` is *not* a log: for `x = TS*p` small it is
approximately linear, so `delta` is an **absolute** difference in proportion; for `x` large it
is logarithmic, so `delta` is a **ratio**. The crossover sits at x ~ 1, i.e. **20 CPM** — right
in the middle of the tested gene universe. The local compression factor `x/(1+x)` is 0.2 at
5 CPM, 0.5 at 20 CPM, 0.83 at 100 CPM.

So a gene's `delta` depends on **that line's own expression level of that gene**. A gene at
500 CPM in K562 and 20 CPM in the target has its donor delta measured on the near-ratio branch
and its true target delta on the half-compressed branch. Transferring `delta` directly carries
the *donor's* compression into the *target's* space, and the error is systematic, gene-
specific, and largest exactly where most of the tested genes live.

## The fix

Fold changes are the transferable quantity; compression is a property of the readout. So
convert out of the donor's space and back into the target's:

```
lfc      = log(expm1(b_ctrl_donor + delta)) - log(expm1(b_ctrl_donor))
delta_T  = log1p(expm1(b_ctrl_target) * exp(lfc)) - b_ctrl_target
```

This uses the target's control profile — the one signal the Challenge gives us per context,
and the one that has so far failed at everything it was tried on (donor weighting, exp 010;
per-gene modulation, exp 017). Here it is not a heuristic weighting, it is the parameter of a
change of variables that is *known* to be needed.

Variants: average donors in `lfc` space then convert once, or convert each donor then average.
Both, with and without whitening, since the two corrections are independent.
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
EPS = 1e-6


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
    m, s = d.decompose()
    gl = [str(g) for g in d.genes]
    # delta of the *perturbation* is m + s; lfc must be computed from the full delta, since
    # the compression is a property of the absolute level, not of the centred part.
    banks[name] = dict(delta=d.delta, m=m, s=s, b=d.b_ctrl,
                       gpos={g: i for i, g in enumerate(gl)}, gset=set(gl),
                       idx={p: i for i, p in enumerate(d.perts)})
vcc = list(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])
common = [g for g in vcc if all(g in b["idx"] for b in banks.values())]
print(f"panel {len(common)}\n")

rng = np.random.default_rng(0)
res = {}
for held in banks:
    others = [n for n in banks if n != held]
    shared = sorted(set.intersection(*[banks[n]["gset"] for n in [held] + others]))
    feat = [g for g in shared if g not in set(common)]
    col = {n: np.fromiter((banks[n]["gpos"][g] for g in feat), np.int64, len(feat))
           for n in [held] + others}
    bT = banks[held]["b"][col[held]].astype(np.float64)
    eT = np.expm1(bT)
    # Truth: the held-out line's own centred effect, unchanged.
    Sh = banks[held]["s"]
    T = np.vstack([Sh[banks[held]["idx"][p]][col[held]] for p in common]).astype(np.float64)
    print(f"fold {held:14} {len(feat)} genes; target median CPM "
          f"{np.median(1e6 * eT / eT.sum()):.1f}", flush=True)

    # --- incumbent: transfer the centred delta directly, in the donor's own log1p space
    sc = {}
    parts_plain, parts_lfc = {}, {}
    for n in others:
        c = col[n]
        Sn = (banks[n]["s"][:, c] / _norm_scale(banks[n]["s"])[c]).astype(np.float64)
        sc[n] = float(np.median(np.linalg.norm(Sn, axis=1))) or 1.0
        parts_plain[n] = np.vstack([Sn[banks[n]["idx"][p]] for p in common])

        # --- lfc re-expression, from the donor's FULL delta at its own basal level
        bD = banks[n]["b"][c].astype(np.float64)
        eD = np.maximum(np.expm1(bD), EPS)
        D_full = np.vstack([banks[n]["delta"][banks[n]["idx"][p]][c] for p in common]
                           ).astype(np.float64)
        lfc = np.log(np.maximum(np.expm1(bD + D_full), EPS)) - np.log(eD)
        dT = np.log1p(eT[None, :] * np.exp(np.clip(lfc, -10, 10))) - bT[None, :]
        dT = dT - dT.mean(0, keepdims=True)          # centre, as `decompose` does
        dT = dT / np.maximum(np.abs(dT).mean(0), np.percentile(np.abs(dT).mean(0), 25))
        parts_lfc[n] = dT

    ref_p = float(np.median([sc[n] for n in others]))
    P_plain = sum((ref_p / sc[n]) * parts_plain[n] for n in others) / len(others)
    scl = {n: float(np.median(np.linalg.norm(parts_lfc[n], axis=1))) or 1.0 for n in others}
    ref_l = float(np.median(list(scl.values())))
    P_lfc = sum((ref_l / scl[n]) * parts_lfc[n] for n in others) / len(others)

    Bl = []
    for n in others:
        c = col[n]
        Sn = (banks[n]["s"][:, c] / _norm_scale(banks[n]["s"])[c]).astype(np.float64)
        if Sn.shape[0] > 6000:
            Sn = Sn[np.sort(rng.choice(Sn.shape[0], 6000, False))]
        Bl.append((ref_p / sc[n]) * Sn)
    W = Whitener.fit(np.vstack(Bl))

    for label, P in (("plain", P_plain), ("lfc", P_lfc), ("half", 0.5 * (
            P_plain / (np.median(np.linalg.norm(P_plain, axis=1)) or 1.0)
            + P_lfc / (np.median(np.linalg.norm(P_lfc, axis=1)) or 1.0)))):
        for g in (0.0, 0.65):
            res.setdefault(f"{label:6} g={g:g}", {})[held] = pds(W.apply(P, g), T)

print(f"\n{'config':16} " + "  ".join(f"{k.replace('orion_','')[:9]:>10}" for k in banks)
      + f" {'macro':>8}")
out = {}
for k, v in sorted(res.items()):
    macro = float(np.mean(list(v.values())))
    out[k] = {**v, "macro": macro}
    print(f"{k:16} " + "  ".join(f"{v[n]:10.4f}" for n in banks) + f" {macro:8.4f}")
best = max(out, key=lambda k: out[k]["macro"])
print(f"\nbest: {best} -> {out[best]['macro']:.4f}")
Path("experiments/026-lfc-space/results.json").write_text(json.dumps(out, indent=2))
