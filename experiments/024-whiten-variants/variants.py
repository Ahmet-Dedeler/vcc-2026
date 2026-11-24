"""024 — three refinements to the whitening operator.

The spectrum of the fitted covariance explains what the operator is actually doing: the top
eigenvalue is 3.90, the median 0.0036, the minimum 1.7e-4, condition number 23,000. The top
10 directions carry 16% of the variance and the top 100 carry 24%. A ridge of 0.01 -- above
the *median* eigenvalue, so it flattens more than half the spectrum -- changed the result by
0.0003. **The gain therefore comes from suppressing the dominant directions, not from
differentiating among the small ones.** That is the hubness story confirmed from the operator
side, and it suggests three specific refinements.

1. **Per-donor whitening.** Each donor line has its *own* dominant programme, reflecting its
   own growth rate and basal state. Averaging first and whitening once removes the pooled
   programme; whitening each donor by its own covariance first would remove each line's
   private one before they are combined. (Whitening every donor by the *same* operator is
   algebraically identical to whitening the average, so only the per-donor-covariance version
   is a new thing.)

2. **Basis size.** 6,000 perturbations per donor was an arbitrary cap; the donors carry
   9,787 / 17,767 / 17,977. More samples give a better-conditioned covariance exactly where
   the operator is most aggressive.

3. **Percentile flooring** instead of an absolute ridge, so the regularisation is scale-free:
   `w_i = max(ev_i, ev_q)^(-gamma/2)`.
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


def apply_floor(W: Whitener, P: np.ndarray, gamma: float, q: float) -> np.ndarray:
    ev = np.maximum(W.evals, np.percentile(W.evals, q)) if q > 0 else W.evals
    w = ev ** (-gamma / 2.0)
    return ((P @ W.evecs) * w) @ W.evecs.T


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
    full = {n: (ref / sc[n]) * banks[n]["s"][:, col[n]].astype(np.float64) for n in others}

    for cap in (2000, 6000, 20000):
        Bl = []
        for n in others:
            M = full[n]
            if M.shape[0] > cap:
                M = M[np.sort(rng.choice(M.shape[0], cap, False))]
            Bl.append(M)
        W = Whitener.fit(np.vstack(Bl))
        print(f"fold {held:14} cap={cap:5} basis {sum(b.shape[0] for b in Bl)}", flush=True)
        for g in (0.55, 0.65, 0.75):
            res.setdefault(f"pooled cap={cap} g={g:g}", {})[held] = pds(W.apply(P, g), T)
        if cap == 6000:
            for q in (10, 25, 50):
                res.setdefault(f"floor q{q} g=0.65", {})[held] = \
                    pds(apply_floor(W, P, 0.65, q), T)

    # per-donor covariance, whiten each donor's panel block with its own operator
    for g in (0.45, 0.65, 0.85):
        acc = np.zeros_like(P)
        for n in others:
            M = full[n]
            if M.shape[0] > 6000:
                M = M[np.sort(rng.choice(M.shape[0], 6000, False))]
            Wn = Whitener.fit(M)
            acc += 0.5 * Wn.apply(parts[n], g)
        res.setdefault(f"per-donor g={g:g}", {})[held] = pds(acc, T)
    print(f"fold {held:14} per-donor done", flush=True)

print(f"\n{'config':26} " + "  ".join(f"{k.replace('orion_','')[:9]:>10}" for k in banks)
      + f" {'macro':>8}")
out = {}
for k, v in sorted(res.items()):
    macro = float(np.mean(list(v.values())))
    out[k] = {**v, "macro": macro}
    print(f"{k:26} " + "  ".join(f"{v[n]:10.4f}" for n in banks) + f" {macro:8.4f}")
best = max(out, key=lambda k: out[k]["macro"])
print(f"\nbest: {best} -> {out[best]['macro']:.4f}")
Path("experiments/024-whiten-variants/results.json").write_text(json.dumps(out, indent=2))
