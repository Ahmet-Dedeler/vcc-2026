"""038 — whiten the genes the operator currently cannot reach.

Whitening is the largest lever in the project (+0.070 macro PDS). It is fitted on genes **every**
donor measures — 7,312 of the 18,533 on the official axis. But the scorer tests 9,923 genes in
context A at the CPM>5 gate, of which **3,131 are off that axis** (Orion-only), and the model
predicts on all of them. So roughly a third of the tested universe is passed through untouched
under our biggest transform.

Experiment 025 fixed the *scaling* consequence (`--whiten-renorm preserve`, +0.006). This asks
the bigger question: does whitening those genes too actually help?

## The test has to be run on the production gene axis

Every LODO fold so far scored on the **donor intersection**, where the situation cannot arise by
construction. Here the folds are built on the **official 18,533-gene axis**, exactly as the
submission is: the feature space is every gene the held-out line shares with the target axis, and
the operator is fitted on either

  `inter`  genes every predictor carries        (what v6/v8 ship — 7.5k genes)
  `full`   genes ANY predictor carries          (~18k genes, so ~98% of the tested universe)

Holding out an Orion line is the informative fold, because K562-gwps carries only 7,679 target
genes while Orion carries 18,106 — so `full` reaches nearly everything and `inter` reaches 41%.

An 18k x 18k eigendecomposition is ~5e13 flops and a few GB; that is the whole reason this was
not tested earlier.
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
GAMMA = [0.3, 0.4, 0.5, 0.6, 0.7]
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
    # column map onto the OFFICIAL axis, which is what production uses
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
print(f"panel {len(common)}; official axis {len(TARGET)} genes")
for n, b in banks.items():
    print(f"  {n:14} covers {len(b['tcols']):,} of the official axis")

rng = np.random.default_rng(0)
res = {}
for held in ("orion_hct116", "orion_hek293t"):        # the informative folds
    others = [n for n in banks if n != held]
    feat = sorted((banks[held]["tset"] - PANEL))
    fidx = {g: i for i, g in enumerate(feat)}
    hmap = {g: j for g, j in zip(np.asarray(TARGET)[banks[held]["tcols"]],
                                 banks[held]["dcols"])}
    T = np.vstack([banks[held]["raw"][banks[held]["idx"][p]][[hmap[g] for g in feat]]
                   for p in common]).astype(np.float64)

    sc, blocks, cols = {}, {}, {}
    for n in others:
        gmap = {g: j for g, j in zip(np.asarray(TARGET)[banks[n]["tcols"]], banks[n]["dcols"])}
        ok = [g for g in feat if g in gmap]
        cols[n] = (np.array([fidx[g] for g in ok]), np.array([gmap[g] for g in ok]))
        sc[n] = float(np.median(np.linalg.norm(banks[n]["s"][:, cols[n][1]], axis=1))) or 1.0
    ref = float(np.median(list(sc.values())))
    P = np.zeros((len(common), len(feat)))
    for n in others:
        ci, di = cols[n]
        P[:, ci] += (1.0 / len(others)) * (ref / sc[n]) * np.vstack(
            [banks[n]["s"][banks[n]["idx"][p]][di] for p in common])

    inter = set(feat)
    for n in others:
        inter &= set(np.asarray(TARGET)[banks[n]["tcols"]])
    axes = {"full": np.array(sorted(set().union(*[set(cols[n][0]) for n in others])))}
    print(f"\nfold {held}: {len(feat):,} scored genes | axes "
          f"{ {k: len(v) for k, v in axes.items()} }", flush=True)
    res.setdefault("none", {})[held] = pds(P, T)

    for aname, ax in axes.items():
        B = np.zeros((0, len(ax)))
        parts = []
        for n in others:
            ci, di = cols[n]
            sub = np.zeros((banks[n]["s"].shape[0], len(feat)), np.float32)
            sub[:, ci] = banks[n]["s"][:, di]
            keep = sub[:, ax]
            if keep.shape[0] > 12000:
                keep = keep[np.sort(rng.choice(keep.shape[0], 12000, False))]
            parts.append((ref / sc[n]) * keep.astype(np.float64))
            del sub
        B = np.vstack(parts); del parts
        print(f"  fitting {aname:6} operator on {B.shape} ...", flush=True)
        W = Whitener.fit(B); del B
        for g in GAMMA:
            Pg = P.copy()
            before = float(np.median(np.linalg.norm(P[:, ax], axis=1)))
            Pg[:, ax] = W.apply(P[:, ax], g, 0.0, renorm=False)
            Pg[:, ax] *= before / (float(np.median(np.linalg.norm(Pg[:, ax], axis=1))) or 1.0)
            res.setdefault(f"{aname} g={g:g}", {})[held] = pds(Pg, T)
            print(f"    {aname} gamma={g}: pds {res[f'{aname} g={g:g}'][held]:.4f}", flush=True)
        del W

print(f"\n{'config':16} " + "  ".join(f"{k.replace('orion_','')[:9]:>10}" for k in
                                      ("orion_hct116", "orion_hek293t")) + f" {'macro':>8}")
out = {}
for k, v in sorted(res.items()):
    macro = float(np.mean(list(v.values())))
    out[k] = {**v, "macro": macro}
    print(f"{k:16} " + "  ".join(f"{v[n]:10.4f}" for n in
                                 ("orion_hct116", "orion_hek293t")) + f" {macro:8.4f}")
best = max(out, key=lambda k: out[k]["macro"])
print(f"\nbest: {best} -> {out[best]['macro']:.4f}")
Path("experiments/038-full-axis/results.json").write_text(json.dumps(out, indent=2))
