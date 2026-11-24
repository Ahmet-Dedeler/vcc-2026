"""025 — what whitening does to the genes it cannot reach.

The whitening operator is fitted on genes **every** donor measures (7,312 of 18,533 on the
production axis). But the model predicts on every gene **any** donor measures, and at the
CPM>20 gate the scorer tests 7,989 genes of which 6,102 are on the operator axis and 1,852
are not — Orion-only genes, ~23% of the tested universe.

Whitening multiplies the on-axis block's norm by ~4.6x, and the builder then rescales the
whole matrix back to its original median row norm. Net effect: **the off-axis genes are
silently down-weighted 4.6x relative to the on-axis ones.** That may be fine — those genes
have a single source, no cross-donor consensus behind them — but it was never measured,
because every LODO fold so far scored on the donor *intersection*, where the situation cannot
arise.

This fold construction reproduces it. The feature axis is every gene the held-out line shares
with **any** predictor; the operator axis is the subset **all** predictors carry. Holding out
Orion-HEK293T with predictors gwps + Orion-HCT116 leaves thousands of scored genes off the
operator axis, exactly as in production.

Three treatments:
  none      no whitening
  global    whiten on-axis, then rescale the whole matrix (what v5 does)
  preserve  whiten on-axis, rescale the on-axis block back to its own prior norm, leaving
            off-axis genes at their original weight
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

rng = np.random.default_rng(0)
res = {}
for held in banks:
    others = [n for n in banks if n != held]
    feat = sorted((banks[held]["gset"] & set.union(*[banks[n]["gset"] for n in others]))
                  - set(common))
    opset = set.intersection(*[banks[n]["gset"] for n in others]) & set(feat)
    fidx = {g: i for i, g in enumerate(feat)}
    opcols = np.array(sorted(fidx[g] for g in opset), dtype=np.int64)
    hcol = np.fromiter((banks[held]["gpos"][g] for g in feat), np.int64, len(feat))
    T = np.vstack([banks[held]["raw"][banks[held]["idx"][p]][hcol]
                   for p in common]).astype(np.float64)

    # Production-shaped prediction: every donor contributes only on the genes it carries,
    # and the 1/n_donors weight is NOT renormalised per gene -- a gene one donor lacks is
    # diluted, exactly as `MultiDonorModel.effect` does it.
    sc, P, Bfull = {}, np.zeros((len(common), len(feat))), {}
    for n in others:
        cols_in_feat = np.array([fidx[g] for g in feat if g in banks[n]["gset"]], np.int64)
        dcols = np.fromiter((banks[n]["gpos"][g] for g in feat if g in banks[n]["gset"]),
                            np.int64, len(cols_in_feat))
        sc[n] = float(np.median(np.linalg.norm(banks[n]["s"][:, dcols], axis=1))) or 1.0
        Bfull[n] = (cols_in_feat, dcols)
    ref = float(np.median(list(sc.values())))
    for n in others:
        cols_in_feat, dcols = Bfull[n]
        block = np.vstack([banks[n]["s"][banks[n]["idx"][p]][dcols] for p in common])
        P[:, cols_in_feat] += (1.0 / len(others)) * (ref / sc[n]) * block

    # Operator basis: all predictor perturbations, on the operator axis only.
    Bl = []
    for n in others:
        cols_in_feat, dcols = Bfull[n]
        keep = np.isin(cols_in_feat, opcols)
        M = banks[n]["s"][:, dcols[keep]]
        if M.shape[0] > 6000:
            M = M[np.sort(rng.choice(M.shape[0], 6000, False))]
        Bl.append((ref / sc[n]) * M.astype(np.float64))
    order = np.argsort(np.concatenate([Bfull[n][0][np.isin(Bfull[n][0], opcols)]
                                       for n in others])[:0]) if False else None
    W = Whitener.fit(np.vstack(Bl))
    off = P.shape[1] - len(opcols)
    print(f"fold {held:14} {P.shape[1]} scored genes, {len(opcols)} on operator axis, "
          f"{off} off ({off/P.shape[1]:.0%})", flush=True)

    res.setdefault("none", {})[held] = pds(P, T)
    for g in (0.45, 0.65, 0.85):
        Pg = P.copy()
        before_all = float(np.median(np.linalg.norm(P, axis=1)))
        before_op = float(np.median(np.linalg.norm(P[:, opcols], axis=1)))
        Pg[:, opcols] = W.apply(P[:, opcols], g, renorm=False)
        after_op = float(np.median(np.linalg.norm(Pg[:, opcols], axis=1))) or 1.0
        glob = Pg * (before_all / (float(np.median(np.linalg.norm(Pg, axis=1))) or 1.0))
        res.setdefault(f"global g={g:g}", {})[held] = pds(glob, T)
        pres = Pg.copy()
        pres[:, opcols] *= before_op / after_op
        res.setdefault(f"preserve g={g:g}", {})[held] = pds(pres, T)

print(f"\n{'config':20} " + "  ".join(f"{k.replace('orion_','')[:9]:>10}" for k in banks)
      + f" {'macro':>8}")
out = {}
for k, v in sorted(res.items()):
    macro = float(np.mean(list(v.values())))
    out[k] = {**v, "macro": macro}
    print(f"{k:20} " + "  ".join(f"{v[n]:10.4f}" for n in banks) + f" {macro:8.4f}")
best = max(out, key=lambda k: out[k]["macro"])
print(f"\nbest: {best} -> {out[best]['macro']:.4f}")
Path("experiments/025-offaxis/results.json").write_text(json.dumps(out, indent=2))
