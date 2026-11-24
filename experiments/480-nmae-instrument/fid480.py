"""480 resume-2: would 434f's iid sign-accuracy reading of the rivals' nmae also show up in their fid?

fid (metrics BRIEF §4): F_p = k / max(n_pred, n_real); n_real = the reference's significant genes (own target excluded);
the call set = our significant genes that the reference can adjudicate (finite, non-zero reference lfc), own target
excluded; k = calls whose predicted and reference lfc share a sign. Perturbations with n_real = n_pred = 0 are dropped;
with n_real = 0 < n_pred, F_p = k / n_pred. fid = mean over the remaining perturbations.

Known answer: 435's H1 harness reports raw fid 0.5089291 for the shipped coarse-MDE arm
(`435/work/h1/metrics-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.json`, computed by cell_eval2's own
dispatch_de_metrics); this script must reproduce it from the same stats (p, q, lfc) without importing cell_eval2.

Experiment (the §5.2 construction): raise our |l|-weighted gate sign accuracy by +Δa by flipping random wrong-signed gate
cells to right, iid (20 draws). The same flipped signs are applied to the predicted lfc; our call set is unchanged
(calls come from the prediction's Wilcoxon, which a sign flip of the arithmetic mean need not move). Read fid and nmae.
Two scopes: flips on gate cells only (the nmae-relevant cells), and flips on every adjudicable called cell (an upper
bound for how much of the accuracy gain lands on our calls).

    uv run python fid480.py      # -> work/fid480.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
keep = ref["keep"].astype(bool)
genes = ref["genes"].astype(str)[keep]
targets = ref["targets"].astype(str)
rq, rl = ref["q"].astype(float), ref["lfc"].astype(float)
own = genes[None, :] == targets[:, None]
st = np.load(ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz",
             allow_pickle=True)
assert (st["targets"].astype(str) == targets).all()
pq, pl_ = st["q"].astype(float), st["lfc"].astype(float)

real_sig = (rq < 0.05) & ~own
adjud = np.isfinite(rl) & (rl != 0)
called = (pq < 0.05) & adjud & ~own


def fid(pred_lfc: np.ndarray) -> float:
    k = (called & (np.sign(pred_lfc) == np.sign(rl))).sum(1)
    n_pred, n_real = called.sum(1), real_sig.sum(1)
    keep_p = ~((n_pred == 0) & (n_real == 0))
    den = np.maximum(n_pred, n_real)[keep_p]
    return float(np.mean(k[keep_p] / den))


gate = real_sig & adjud & np.isfinite(rl)
scored = gate.sum(1) >= 10
gate &= scored[:, None]
Tz = np.where(gate, rl, 0.0)


def nmae(pred: np.ndarray) -> float:
    P = np.where(np.isfinite(pred), pred, 0.0)
    return float(np.mean((np.abs(P - Tz) * gate).sum(1)[scored] / np.abs(Tz).sum(1)[scored]))


f0 = fid(pl_)
print(f"known answer: fid {f0:.7f} (harness 0.5089291)   nmae {nmae(pl_):.5f} (0.94270)")
if abs(f0 - 0.5089291) > 2e-6:
    print("KNOWN ANSWER FAILED")
    sys.exit(1)
n_pred, n_real = called.sum(1), real_sig.sum(1)
yield_p = np.minimum(1, n_pred / np.maximum(n_real, 1))
prec = (called & (np.sign(pl_) == np.sign(rl))).sum(1) / np.maximum(n_pred, 1)
share_called = float((np.abs(Tz) * (gate & called)).sum() / np.abs(Tz).sum())
print(f"median n_pred {np.median(n_pred):.0f}, median n_real {np.median(n_real):.0f}, mean yield {yield_p.mean():.3f}, "
      f"mean precision (targets with calls) {prec[n_pred > 0].mean():.3f}; called share of gate |l| {share_called:.3f}")

w = np.abs(Tz)
right = (np.sign(pl_) == np.sign(rl)) & gate
acc0 = float((w * right).sum() / w.sum())
out = dict(ka_fid=f0, n_pred_median=float(np.median(n_pred)), n_real_median=float(np.median(n_real)),
           yield_mean=float(yield_p.mean()), gate_acc=acc0, called_share_gate_mass=share_called, rows=[])
print(f"gate |l|-weighted sign accuracy {acc0:.4f}")
rng = np.random.default_rng(0)
for scope in ["gate", "gate+calls"]:
    pool = gate & ~right & (pl_ != 0)
    if scope == "gate+calls":
        pool_calls = called & ~gate & (np.sign(pl_) != np.sign(rl)) & (pl_ != 0)
    for da in [0.02, 0.05, 0.08]:
        fs, ns = [], []
        for _ in range(20):
            idx = np.flatnonzero(pool.ravel())
            rng.shuffle(idx)
            wi = w.ravel()[idx]
            m = int(np.searchsorted(np.cumsum(wi), da * w.sum()))       # flip until +da of |l| mass is right
            P = pl_.copy().ravel()
            P[idx[:m]] = -P[idx[:m]]
            if scope == "gate+calls":
                # the same flip RATE (share of wrong cells turned right) on called, adjudicable, non-gate cells
                rate = m / max(idx.size, 1)
                jdx = np.flatnonzero(pool_calls.ravel())
                rng.shuffle(jdx)
                P[jdx[:int(round(rate * jdx.size))]] *= -1
            P = P.reshape(pl_.shape)
            fs.append(fid(P))
            ns.append(nmae(P))
        r = dict(scope=scope, da=da, d_fid=float(np.mean(fs) - f0), d_fid_sd=float(np.std(fs)),
                 d_nmae=float(np.mean(ns) - nmae(pl_)))
        out["rows"].append(r)
        print(f"{scope:10s} +{da:.2f} acc: Δfid {r['d_fid']:+.4f} (sd {r['d_fid_sd']:.4f})  Δnmae {r['d_nmae']:+.4f}")
json.dump(out, open(HERE / "work/fid480.json", "w"), indent=1)
