"""039 — genes only some donors cover are silently shrunk. Should they be?

`MultiDonorModel.effect` divides by the total weight over donors carrying the *perturbation*,
then adds each donor's contribution only on the genes that donor measures. So for a gene that
K562-gwps does not carry — and gwps carries 7,679 of the official axis where each Orion line
carries 18,111 — the prediction is

    (1/3)*0 + (1/3)*hct116 + (1/3)*hek293t  =  2/3 of the Orion average

**Two thirds of the gene axis is emitted at 2/3 magnitude**, not because the estimate is weaker
but because a donor that never measured the gene is being counted as a vote of zero.

This did not matter much when the whitening operator only covered the intersection anyway.
Experiment 038 changed that: the operator now spans 17,942 genes, so those columns carry real
weight in the direction.

Whether to fix it is genuinely ambiguous, which is why it needs measuring rather than arguing:

  as-is       under-covered genes are shrunk toward zero. That is accidental, but it is also
              exactly what a confidence weighting would do -- a 2-donor mean is noisier than a
              3-donor mean, and shrinking noisy estimates is what every other part of this
              pipeline tries to do deliberately.
  per-gene    divide each gene by the donors that actually cover it, so every gene gets an
              unbiased mean of whatever measured it.

Run on the production gene axis, as experiment 038 is, with the union operator.
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
TARGET = list(pd.read_csv("data/challenge/gene_names.csv")["gene_name"])
TIDX = {g: i for i, g in enumerate(TARGET)}
PANEL = set(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])


def pds(P, T):
    Pn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12)
    Tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    D = 1.0 - Pn @ Tn.T
    n = D.shape[1]
    sd = np.diag(D)[:, None]
    return float((1.0 - ((D < sd).sum(1) + 0.5 * ((D == sd).sum(1) - 1)) / (n - 1)).mean())


banks = {}
for name, path in DONORS.items():
    d = Donor.load(path)
    _, s = d.decompose()
    gl = [str(g) for g in d.genes]
    tc, dc = [], []
    for j, g in enumerate(gl):
        t = TIDX.get(g)
        if t is not None:
            tc.append(t); dc.append(j)
    banks[name] = dict(s=(s / _norm_scale(s)).astype(np.float32), raw=s,
                       tcols=np.asarray(tc), dcols=np.asarray(dc),
                       tset=set(np.asarray(TARGET)[tc]),
                       idx={p: i for i, p in enumerate(d.perts)})
common = [g for g in PANEL if all(g in b["idx"] for b in banks.values())]
rng = np.random.default_rng(0)
res = {}
for held in ("orion_hct116", "orion_hek293t"):
    others = [n for n in banks if n != held]
    feat = sorted(banks[held]["tset"] - PANEL)
    fidx = {g: i for i, g in enumerate(feat)}
    hmap = {g: j for g, j in zip(np.asarray(TARGET)[banks[held]["tcols"]], banks[held]["dcols"])}
    T = np.vstack([banks[held]["raw"][banks[held]["idx"][p]][[hmap[g] for g in feat]]
                   for p in common]).astype(np.float64)
    sc, cols = {}, {}
    for n in others:
        gmap = {g: j for g, j in zip(np.asarray(TARGET)[banks[n]["tcols"]], banks[n]["dcols"])}
        ok = [g for g in feat if g in gmap]
        cols[n] = (np.array([fidx[g] for g in ok]), np.array([gmap[g] for g in ok]))
        sc[n] = float(np.median(np.linalg.norm(banks[n]["s"][:, cols[n][1]], axis=1))) or 1.0
    ref = float(np.median(list(sc.values())))
    P = np.zeros((len(common), len(feat)))
    cover = np.zeros(len(feat))
    for n in others:
        ci, di = cols[n]
        P[:, ci] += (ref / sc[n]) * np.vstack(
            [banks[n]["s"][banks[n]["idx"][p]][di] for p in common])
        cover[ci] += 1.0
    n_full = int((cover == len(others)).sum())
    print(f"fold {held}: {len(feat):,} genes | covered by all {len(others)} predictors "
          f"{n_full:,} ({100*n_full/len(feat):.0f}%)", flush=True)

    variants = {"as-is (divide by n_donors)": P / len(others),
                "per-gene (divide by coverage)": P / np.maximum(cover, 1)[None, :]}
    ax = np.array(sorted(set().union(*[set(cols[n][0]) for n in others])))
    parts = []
    for n in others:
        ci, di = cols[n]
        sub = np.zeros((banks[n]["s"].shape[0], len(feat)), np.float32)
        sub[:, ci] = banks[n]["s"][:, di]
        parts.append((ref / sc[n]) * sub[:, ax].astype(np.float64)); del sub
    W = Whitener.fit(np.vstack(parts)); del parts
    for vname, Pv in variants.items():
        for g in (0.0, 0.65):
            Pg = Pv.copy()
            if g:
                b = float(np.median(np.linalg.norm(Pv[:, ax], axis=1)))
                Pg[:, ax] = W.apply(Pv[:, ax], g, 0.0, renorm=False)
                Pg[:, ax] *= b / (float(np.median(np.linalg.norm(Pg[:, ax], axis=1))) or 1.0)
            res.setdefault(f"{vname} g={g:g}", {})[held] = pds(Pg, T)
            print(f"  {vname:32} gamma={g}: {res[f'{vname} g={g:g}'][held]:.4f}", flush=True)
    del W

print(f"\n{'config':40} {'hct116':>9} {'hek293t':>9} {'macro':>9}")
out = {}
for k, v in sorted(res.items()):
    m = float(np.mean(list(v.values())))
    out[k] = {**v, "macro": m}
    print(f"{k:40} {v['orion_hct116']:9.4f} {v['orion_hek293t']:9.4f} {m:9.4f}")
best = max(out, key=lambda k: out[k]["macro"])
print(f"\nbest: {best} -> {out[best]['macro']:.4f}")
Path("experiments/039-coverage-weight/results.json").write_text(json.dumps(out, indent=2))
