"""021 — whiten by the NOISE covariance, not the signal covariance.

Experiment 017 whitens with `C_signal = B^T B / n`, the covariance of donor effect vectors.
That is the standard "flatten the dominant programmes" argument, and it works. But it is not
what the estimation theory actually asks for.

Write a prediction as `P = T + e`: the true target effect plus an error made of measurement
noise and line-specific biology. For a matched-filter / retrieval problem the optimal
similarity is Mahalanobis in the **error** covariance, not the signal covariance:

    similarity(P, T) ~ P^T Sigma_e^{-1} T

Whitening by `C_signal` is only equivalent to that when the error is proportional to the
signal. The measured optimum sitting at gamma = 0.65 rather than 1.0 is consistent with the
truth being somewhere between "error is isotropic" (gamma = 0) and "error is proportional to
signal" (gamma = 1) -- which is a hint that the right object has not been used yet.

**Sigma_e is directly estimable.** Two donors measuring the same perturbation differ by
exactly the thing we are calling error: sampling noise plus line-specific response. So

    Sigma_e ~ cov over perturbations of (d_A - d_B) / 2

computed on every perturbation the two predictor lines share -- thousands of them, and no
target data anywhere. Whether the signal basis or the noise basis wins is an empirical
question about how the error is structured, and it has never been asked here.

Three bases compared, all at matched gamma grids:
  signal  C = B^T B / n over all donor effects            (experiment 017's operator)
  noise   C = D^T D / n over donor-difference vectors
  ratio   signal-whiten then noise-colour, i.e. the generalised eigenproblem direction
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
GAMMAS = [0.0, 0.3, 0.45, 0.65, 0.85, 1.1]


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
print(f"panel {len(common)}")

rng = np.random.default_rng(0)
results = {}
for held in banks:
    others = [n for n in banks if n != held]
    shared_g = sorted(set.intersection(*[banks[n]["gset"] for n in [held] + others]))
    feat = [g for g in shared_g if g not in set(common)]
    col = {n: np.fromiter((banks[n]["gpos"][g] for g in feat), np.int64, len(feat))
           for n in [held] + others}
    T = np.vstack([banks[held]["raw"][banks[held]["idx"][p]][col[held]]
                   for p in common]).astype(np.float64)
    sc = {n: float(np.median(np.linalg.norm(banks[n]["s"][:, col[n]], axis=1))) or 1.0
          for n in others}
    ref = float(np.median(list(sc.values())))
    parts = {n: (ref / sc[n]) * np.vstack(
        [banks[n]["s"][banks[n]["idx"][p]][col[n]] for p in common]).astype(np.float64)
        for n in others}
    P = 0.5 * (parts[others[0]] + parts[others[1]])

    # Signal basis: every perturbation each predictor measured.
    Bl = []
    for n in others:
        M = banks[n]["s"][:, col[n]]
        if M.shape[0] > 6000:
            M = M[np.sort(rng.choice(M.shape[0], 6000, False))]
        Bl.append((ref / sc[n]) * M.astype(np.float64))
    # Noise basis: half-differences on the perturbations both predictors carry.
    both = [p for p in set(banks[others[0]]["idx"]) & set(banks[others[1]]["idx"])]
    if len(both) > 6000:
        both = list(np.asarray(both)[np.sort(rng.choice(len(both), 6000, False))])
    A = (ref / sc[others[0]]) * np.vstack(
        [banks[others[0]]["s"][banks[others[0]]["idx"][p]][col[others[0]]] for p in both])
    Bm = (ref / sc[others[1]]) * np.vstack(
        [banks[others[1]]["s"][banks[others[1]]["idx"][p]][col[others[1]]] for p in both])
    Dm = 0.5 * (A.astype(np.float64) - Bm.astype(np.float64))
    print(f"fold {held:14} signal basis {sum(b.shape[0] for b in Bl)}, "
          f"noise basis {Dm.shape[0]} shared perturbations", flush=True)

    W_sig = Whitener.fit(np.vstack(Bl))
    W_noi = Whitener.fit(Dm)
    for label, W in (("signal", W_sig), ("noise", W_noi)):
        for g in GAMMAS:
            results.setdefault(f"{label} gamma={g:g}", {})[held] = pds(W.apply(P, g), T)
    # ratio: noise-whiten, then re-colour by signal -- the generalised-eigenproblem direction
    for g in GAMMAS[1:]:
        Pw = W_sig.apply(W_noi.apply(P, g), -g)
        results.setdefault(f"ratio gamma={g:g}", {})[held] = pds(Pw, T)

print(f"\n{'config':22} " + "  ".join(f"{k.replace('orion_','')[:9]:>10}" for k in banks)
      + f" {'macro':>8}")
out = {}
for k, v in results.items():
    macro = float(np.mean(list(v.values())))
    out[k] = {**v, "macro": macro}
    print(f"{k:22} " + "  ".join(f"{v[n]:10.4f}" for n in banks) + f" {macro:8.4f}")
best = max(out, key=lambda k: out[k]["macro"])
print(f"\nbest: {best} -> {out[best]['macro']:.4f}")
Path("experiments/021-noise-whiten/results.json").write_text(json.dumps(out, indent=2))
