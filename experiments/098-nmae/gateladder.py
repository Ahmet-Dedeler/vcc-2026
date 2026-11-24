"""THE transfer question: the challenge scores ~4,757 reference-significant genes per row and
our harnesses score 89 (orion) / 1,119 (h1).  Does the arm ranking depend on gate SIZE?

a13 is the test case with a known answer.  Harness DE gates said its `nmae` moved by
+0.0048/+0.0180 raw (h1 said BETTER on the builder's own form, -0.0064); the board delivered
**+0.0744**.  If the miss is gate composition, the penalty must grow monotonically as the gate
is relaxed toward the challenge's depth.

The gate is relaxed by raising the BH threshold, which is a PROXY for a more powerful reference
(more, weaker genes), not the same thing.  Labelled derived, not measured.
"""
from __future__ import annotations
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/098-nmae"); sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
import nmaekit as K
from linear import donor_keys
from scipy.stats import rankdata

ALPHAS = [0.05, 0.15, 0.35, 0.6, 0.85, 1.01]


def gates(which, alpha):
    z = np.load(f"{ROOT}/experiments/098-nmae/work/refde_{which}.npz", allow_pickle=True)
    q = z["q"]; lf = z["lfc"].astype(np.float64); perts = z["perts"].astype(str)
    out = {}
    for i, p in enumerate(perts):
        s = (q[i] < alpha) & np.isfinite(lf[i]) & (lf[i] != 0)
        out[p] = np.flatnonzero(s)
    return out, lf, perts


def desc_midrank(x):
    n = x.size; o = np.argsort(-x, kind="stable"); od = x[o]
    e = np.flatnonzero(np.concatenate(([True], od[1:] != od[:-1], [True])))
    st, en = e[:-1], e[1:]
    out = np.empty(n); out[o] = np.repeat((st + en - 1) / 2.0, en - st)
    return out


for W_ in ("h1", "orion"):
    H = K.Harness(W_); eff0, _ = K.a12_eff(H)
    ag, _ = donor_keys(H, W_)
    g = H.gate
    # a13 as BUILT: alpha_g = exp(-gamma * descending midrank(agreement) / 2000), gamma = 0.25,
    # times the per-perturbation zr**1 factor, at its own re-solved cliff amplitude.
    A13 = np.ones_like(eff0)
    for i in range(H.P):
        A13[i] = 1.0 if np.unique(ag[i]).size < 2 else np.exp(-0.25 * desc_midrank(ag[i]) / 2000.0)
    from cell_eval2.metrics.delta import _across_pert_budget
    B13 = float(_across_pert_budget(A13 * H.dh0, H.cols, H.w))
    s13 = float(np.sqrt(H.claim / B13))
    eff13 = s13 * A13 * H.dh0
    L0 = H.lfc_hat(eff0[H.rows]); L13 = H.lfc_hat(eff13[H.rows])
    Lc = {c: H.lfc_hat(eff0[H.rows]) * 0 for c in ()}
    rpos = {p: i for i, p in enumerate(np.load(f"{ROOT}/experiments/098-nmae/work/refde_{W_}.npz",
                                               allow_pickle=True)["perts"].astype(str))}
    print(f"\n===== {W_}  (a13 cliff s = {s13:.4f}, mean shrink weight {A13.mean():.4f}) =====")
    print(f"  {'BH alpha':>9}{'mean |S_p|':>12}{'rows':>6}{'nmae a12':>10}{'nmae a13':>10}"
          f"{'Δa13':>9}{'c* (lfc)':>10}{'nmae@c*':>10}{'Δc*':>9}{'sign acc':>10}")
    for al in ALPHAS:
        gt, lf, perts = gates(W_, al)
        rows, rr = [], []
        for k, (i, j) in enumerate(zip(H.rows, H.rrows)):
            gs = gt[perts[j]]
            if gs.size >= 10: rows.append(k); rr.append((j, gs))
        if not rows: continue
        n0 = n13 = 0.0; sz = []; acc = []; hs, rs_ = [], []
        for k, (j, gs) in zip(rows, rr):
            r = lf[j][gs]; den = np.abs(r).sum(); sz.append(gs.size)
            n0 += np.abs(L0[k][gs] - r).sum() / den
            n13 += np.abs(L13[k][gs] - r).sum() / den
            hs.append(L0[k][gs]); rs_.append(r)
            acc.append(np.average(np.sign(L0[k][gs]) == np.sign(r),
                                  weights=np.abs(L0[k][gs]) + 1e-12))
        n0 /= len(rows); n13 /= len(rows)
        cs = np.arange(0.2, 6.01, 0.02)
        f = np.array([np.mean([np.abs(c * a - b).sum() / np.abs(b).sum()
                               for a, b in zip(hs, rs_)]) for c in cs])
        bi = int(f.argmin())
        print(f"  {al:>9.2f}{np.mean(sz):>12.0f}{len(rows):>6}{n0:>10.4f}{n13:>10.4f}"
              f"{n13-n0:>+9.4f}{cs[bi]:>10.2f}{f[bi]:>10.4f}{f[bi]-n0:>+9.4f}"
              f"{np.mean(acc):>10.4f}")
