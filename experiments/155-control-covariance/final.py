"""155 — the consolidated table that RESULTS.md reports.  One frame, one run, no re-derivation.

Both sign conventions are scored explicitly rather than asserted from `acc(+b) = 1 - acc(-b)`:
the identity holds only where no coefficient is exactly zero, and the argument the experiment
turns on should not rest on that.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bench import HDR, ROW, load_frame, score  # noqa: E402

s, inc, gdst, pok = load_frame()


def frame(M):
    P = np.zeros_like(s.pred_lfc)
    P[:, gdst] = M[:, gdst]
    return P


def dbl(M):
    M = M - M.mean(0, keepdims=True)
    return M - M.mean(1, keepdims=True)


def match(P):
    a = np.linalg.norm(inc[:, gdst], axis=1)
    b = np.maximum(np.linalg.norm(P[:, gdst], axis=1), 1e-30)
    return P * (a / b)[:, None]


def main() -> None:
    print(f"frame: {pok.sum()} perturbations x {gdst.size} genes; reference DE pool "
          f"{sum(int(s.real_sig[p].sum()) for p in range(len(s.perts)))} gene-calls, "
          f"59.7% down\n")
    print(HDR)
    base = score(s, inc)
    print(ROW.format(name="INCUMBENT 104 baseline g=0.65", **base))

    rng = np.random.default_rng(0)
    ad = frame(-np.ones((len(s.perts), s.genes.size)))
    ad[:, gdst] *= rng.random(gdst.size)
    print(ROW.format(name="null: all-down constant", **score(s, ad)))
    print()

    B = {k: np.load(os.path.join(HERE, "work", f"coef_k{k}.npz"),
                    allow_pickle=True)["beta"].astype(np.float64) for k in (0, 5, 20, 50)}
    for k in (0, 5, 20, 50):
        for sgn, tag in ((-1, "-beta  (knockdown = T down)"), (+1, "+beta  (sign flipped)")):
            P = match(frame(sgn * dbl(B[k])))
            print(ROW.format(name=f"k={k:<2} {tag}", **score(s, P)))
        print()

    print("BLENDS  (1-f)*incumbent + f*candidate, candidate L2-matched per row")
    print(HDR + "   d(sign acc)  d(head)   d(reach)   d(pds)")
    for k, sgn, lab in ((0, -1, "k=0 -beta"), (0, +1, "k=0 +beta"), (20, -1, "k=20 -beta")):
        C = match(frame(sgn * dbl(B[k])))
        for f in (0.02, 0.05, 0.10, 0.20):
            r = score(s, (1 - f) * inc + f * C)
            print(ROW.format(name=f"{lab}  f={f:.2f}", **r)
                  + f"   {r['sign_acc']-base['sign_acc']:+11.4f} "
                    f"{r['head_purity']-base['head_purity']:+8.4f} "
                    f"{r['reach']-base['reach']:+10.5f} {r['pds']-base['pds']:+8.4f}")
        print()


if __name__ == "__main__":
    main()
