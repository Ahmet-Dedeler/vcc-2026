"""192 Phase 4 — the price, swept over the retrodicted band of the board's `n_real` level, and
the SECOND, independent estimate of that level from a15's own board `jac`.

Two known-answer rows the board supplies for free and that nothing in this repo has used
together before:
    a15 board raw `fid` 0.503536  at challenge n_pred 4,474
    a15 board raw `jac` 0.030080  at the same point
`cover192.py` fits `f` from the first (key-free, through the coverage term).  This script fits
it again from the second (through the union term of `jac`) and reports whether the two agree.
Only then is the arm priced.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
E185 = os.path.join(REPO, "experiments/185-rank-controlled-emission")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, HERE)
from debench import Surface                                            # noqa: E402
from frontier192 import A15, INV, PowerSurface, order_key, topk_call   # noqa: E402
from keys192 import gene_stats                                         # noqa: E402

SHRINK = 0.4277
BOARD_A15 = dict(fid=0.503536, jac=0.030080, reach=0.181173)
NPRED_A15 = 4474
FS = [0.30, 0.40, 0.50, 0.65, 0.80, 1.00]
KS = [1200, 1600, 2200, 2800, 3600, 4474, 6000]


def main() -> None:
    E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
    pos = {str(g): i for i, g in enumerate(E["genes"])}
    tpos = {str(t): i for i, t in enumerate(E["targets"])}
    s = Surface(A15)
    z15 = np.load(A15, allow_pickle=True)
    gi = np.array([pos[g] for g in s.genes])
    ti = np.array([tpos[t] for t in s.perts])
    cpm = E["cpm"][gi].astype(np.float64)
    absH = np.abs(E["H"][np.ix_(ti, gi)].astype(np.float64))
    T = absH.shape[0]
    inc = -np.log(np.maximum(z15["pred_p"], 1e-300))
    gs = gene_stats(os.path.join(E185, "work/ctrl_cs.npy"), os.path.join(WORK, "genestats.npz"))
    sig = s.real_sig & np.isfinite(s.real_lfc)

    from sklearn.linear_model import LogisticRegression
    X = np.stack([np.broadcast_to(np.log(cpm + 1)[None, :], absH.shape),
                  np.broadcast_to(np.log(gs["nz"][:] + 1e-3)[None, :], absH.shape),
                  np.broadcast_to(np.log(gs["cv2"][:] + 1e-3)[None, :], absH.shape)], -1)
    Xz = (X - X.reshape(-1, 3).mean(0)) / X.reshape(-1, 3).std(0)
    rngf = np.random.default_rng(7)
    fitk = np.zeros(absH.shape)
    for p in range(T):
        tr = np.setdiff1d(np.arange(T), [p])
        Xt = Xz[tr].reshape(-1, 3)
        yt = sig[tr].reshape(-1)
        idx = rngf.choice(Xt.shape[0], size=min(300_000, Xt.shape[0]), replace=False)
        fitk[p] = LogisticRegression(max_iter=400).fit(Xt[idx], yt[idx]).decision_function(Xz[p])

    print("=" * 100)
    print("A.  SECOND ESTIMATE OF THE BOARD'S n_real LEVEL, from a15's own board `jac` 0.030080")
    print(f"    (modelled call set = top-{NPRED_A15} by a15's own emitted p-order)")
    print(f"  {'f':>6}{'median n_real':>15}{'model jac':>11}{'board jac':>11}{'ratio':>8}"
          f"{'model fid':>11}{'board fid':>11}{'ratio':>8}")
    base, out = {}, {}
    for f in FS:
        sf = PowerSurface(A15, f)
        nr = np.array([int(sf.real_sig[p].sum()) for p in range(T)])
        call = topk_call(inc, np.full(T, NPRED_A15))
        v = dict(fid=sf.fid(call=call), jac=sf.jac(call=call),
                 reach=sf.reach(key=order_key(call, absH)))
        base[f] = (sf, v)
        print(f"  {f:>6}{np.median(nr):>15.0f}{v['jac']:>11.5f}{BOARD_A15['jac']:>11.5f}"
              f"{v['jac'] / BOARD_A15['jac']:>8.2f}{v['fid']:>11.5f}{BOARD_A15['fid']:>11.5f}"
              f"{v['fid'] / BOARD_A15['fid']:>8.3f}")
    print("    (the h1-test surface is not the challenge -- rule 9.  A ratio near 1 means the "
          "LEVEL transfers;\n     the priced deltas below are read at the f where it does.)")

    print(f"\n{'=' * 100}")
    print("B.  THE PRICED ARM — call set = top-K by the FIT (cpm, nz, cv2) key, LOTO, realised "
          "through 191's emitter.")
    print("    `pds` and `nmae` are PINNED (191 falsifier 1: 900/900 delivered per-gene panel "
          "totals\n     digest-identical; nmae +0.00011 raw).  `mse` is PINNED (185 §6: `claim` "
          "sits at the\n     reference cap 0.24653 for every rc arm).  So SUM(3) below IS the "
          "whole ledger,\n     before rule 13's toll.")
    for f in FS:
        sf, a0 = base[f]
        print(f"\n  f = {f}   a15 baseline: fid {a0['fid']:.5f} jac {a0['jac']:.5f} "
              f"reach {a0['reach']:.5f}")
        print(f"    {'K':>7}{'fid':>9}{'jac':>9}{'reach':>8}{'d fid':>10}{'d jac':>10}"
              f"{'d reach':>10}{'SUM(3)':>10}{'with toll':>11}")
        for K in KS:
            call = topk_call(fitk, np.full(T, K))
            v = dict(fid=sf.fid(call=call), jac=sf.jac(call=call),
                     reach=sf.reach(key=order_key(call, absH)))
            d = {m: (v[m] - a0[m]) * INV[m] * (SHRINK if m == "reach" else 1.0) / 6 for m in v}
            tot = sum(d.values())
            out[f"{f}|{K}"] = dict(raw=v, overall=d, sum3=tot, toll=tot - 0.0146 / 6)
            print(f"    {K:>7}{v['fid']:>9.5f}{v['jac']:>9.5f}{v['reach']:>8.4f}"
                  f"{d['fid']:>+10.5f}{d['jac']:>+10.5f}{d['reach']:>+10.5f}{tot:>+10.5f}"
                  f"{tot - 0.0146 / 6:>+11.5f}")

    with open(os.path.join(WORK, "price192.json"), "w") as fh:
        json.dump(dict(base={str(k): v[1] for k, v in base.items()}, rows=out), fh, indent=1,
                  default=float)
    print(f"\nwrote {WORK}/price192.json")


if __name__ == "__main__":
    main()
