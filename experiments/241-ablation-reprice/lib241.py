"""241 — shared instrument code: the four DE members on h1-test, in overall board units.

Written independently of `six200.py` (the constants are re-read from AGENTS.md / 157 / 169,
not imported from it) so that reproducing 200's published DE legs is a real retrodiction of
200's arithmetic and not a call into it.  The one thing that IS shared is
`debench.Surface` — that is the production instrument every board-linked `reach` number in
`docs/BOARD-LEDGER.md` is on, and replacing it would change what is being measured.

Conventions, all inherited from 169 §4 / 200 §"six200.py" verbatim:
  pool    the reference's own BH<0.05 set (significance-selected).  `n_real_target` is NOT
          used: it deletes the small-`n_real` tercile that carries 98% of our reach.
  reach   `Surface.reach()` × K = (0.1890-0.1811)/(0.2029-0.1828) = 0.3930, the h1-test→board
          transfer measured on 157 §1's own board pair, / span / 6.
  fid     precision `mean_p(k/n_pred)`; the yield cap is saturated on every harness.
  jac     rebuilt on the BOARD's counts from measured recall, n_real 240, n_pred 4424.
  nmae    `Surface.nmae()`, span 0.600, sign flipped (lower raw is better).
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, os.path.join(ROOT, "experiments/153-de-bench"))
from debench import ANCHOR, Surface                                    # noqa: E402

ALPHA = 0.05
K_REACH = (0.1890 - 0.1811) / (0.2029 - 0.1828)     # 157 §1 -> 0.39303...
SPAN_R = ANCHOR["reach"][1] - ANCHOR["reach"][0]
SPAN_F = ANCHOR["fid"][1] - ANCHOR["fid"][0]
SPAN_J = ANCHOR["jac"][1] - ANCHOR["jac"][0]
SPAN_N = 1.0013 - 0.400
BOARD_NR, BOARD_NP = 240.0, 4424.0

# rule 13's composed slot gate
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089


def precision(s: Surface) -> float:
    out = []
    for p in range(len(s.perts)):
        sel = np.flatnonzero(s.pred_sig[p] & s.adj[p])
        if sel.size == 0:
            continue
        out.append(float((np.sign(s.pred_lfc[p, sel]) ==
                          np.sign(s.real_lfc[p, sel])).sum()) / sel.size)
    return float(np.mean(out))


def recall_of(s: Surface) -> float:
    out = []
    for p in range(len(s.perts)):
        nr = int(s.real_sig[p].sum())
        if nr == 0:
            continue
        out.append(int((s.real_sig[p] & s.pred_sig[p]).sum()) / nr)
    return float(np.mean(out))


def row_for(path: str) -> dict:
    s = Surface(path)
    nr, npd = s.counts()
    return dict(n_pred=float(np.median(npd)), n_real=float(np.median(nr)),
                precision=precision(s), recall=recall_of(s),
                reach=float(s.reach()), nmae=float(s.nmae()),
                fid_raw=float(s.fid()), jac_raw=float(s.jac()),
                sign=float(s.sign_accuracy()))


def nmae_parts(s: Surface) -> dict:
    """188's exact `nmae` identity, per perturbation then averaged (BRIEF 7's rule).

        nmae = 1 + U/D - 2 M/D
        D = sum |truth|,  U = sum |pred|,  M = sum over SIGN-MATCHED genes of
        min(|pred|, |truth|)

    `M` can only be collected by getting a SIGN right; a magnitude rebalance moves `U`.
    `e = dM/dU` is therefore the efficiency of an arm, break-even 0.5 (§8.1d).
    """
    P, T = s.pred_lfc, s.real_lfc
    Us, Ms, Ns = [], [], []
    for p in range(len(s.perts)):
        idx = np.flatnonzero(s.real_sig[p] & np.isfinite(T[p]))
        if idx.size < 10:
            continue
        d = np.abs(T[p, idx]).sum()
        if not np.isfinite(d) or d == 0:
            continue
        h = np.where(np.isfinite(P[p, idx]), P[p, idx], 0.0)
        t = T[p, idx]
        m = np.sign(h) == np.sign(t)
        Us.append(np.abs(h).sum() / d)
        Ms.append(np.minimum(np.abs(h[m]), np.abs(t[m])).sum() / d)
        Ns.append(np.abs(h - t).sum() / d)
    U, M, N = float(np.mean(Us)), float(np.mean(Ms)), float(np.mean(Ns))
    return dict(U_over_D=U, M_over_D=M, nmae=N, identity=1.0 + U - 2.0 * M,
                n_scored=len(Us))


def efficiency(arm: dict, base: dict) -> dict:
    """`e = dM/dU` between two `nmae_parts` dicts."""
    dU = arm["U_over_D"] - base["U_over_D"]
    dM = arm["M_over_D"] - base["M_over_D"]
    return dict(dU=dU, dM=dM, e=(dM / dU) if dU != 0 else float("nan"),
                d_nmae=arm["nmae"] - base["nmae"])


def jac_board(recall: float, n_pred: float = BOARD_NP) -> float:
    hit = recall * BOARD_NR
    return hit / (n_pred + BOARD_NR - hit)


def de_overall(row: dict, base: dict) -> dict:
    """One arm's four DE members against its base, in overall board units."""
    dr = row["reach"] - base["reach"]
    b_np = BOARD_NP * (base["n_pred"] / row["n_pred"]) if row.get("n_pred") else BOARD_NP
    return dict(
        d_reach_raw=dr,
        reach=dr * K_REACH / SPAN_R / 6,
        fid=(row["precision"] - base["precision"]) / SPAN_F / 6,
        jac=(jac_board(row["recall"]) - jac_board(base["recall"])) / SPAN_J / 6,
        jac_np=(jac_board(row["recall"], BOARD_NP)
                - jac_board(base["recall"], b_np)) / SPAN_J / 6,
        nmae=-(row["nmae"] - base["nmae"]) / SPAN_N / 6,
        d_nmae_raw=row["nmae"] - base["nmae"])


def paired_de(rows: dict, arm_names, base_names) -> dict:
    """Seed-paired mean of the four members plus the paired error bar on raw `reach`."""
    n = min(len(arm_names), len(base_names))
    per = [de_overall(rows[arm_names[i]], rows[base_names[i]]) for i in range(n)]
    out = {}
    for k in ("reach", "fid", "jac", "jac_np", "nmae", "d_reach_raw", "d_nmae_raw"):
        v = np.asarray([p[k] for p in per])
        out[k] = float(v.mean())
        if k == "d_reach_raw":
            sem = float(v.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
            out["reach_sem_raw"] = sem
            out["reach_t"] = float(v.mean() / sem) if n > 1 and sem > 0 else float("nan")
    out["n_seeds"] = n
    return out


def seed_names(stem: str, n: int):
    return [stem if i == 0 else f"{stem}_s{i + 1}" for i in range(n)]
