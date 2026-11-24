"""155 — is the `+beta` blend gain target-specific, or is it the pool's marginal?

`final.py` shows `+beta` at k=0 standalone reaching 0.1747 (incumbent 0.0904) and a 5%
blend improving sign accuracy, head purity, reach AND pds at once.  That is exactly the
shape of a result this repo has been burned by before, and it has a boring candidate
explanation: `+beta`'s heaviest coefficients sit on the high-variance genes, those genes are
predominantly DOWN in the reference pool (59.7% overall), and `+beta` predicts 55.4% down.
None of that requires the coefficient row to belong to the target it is scored against.

The test: deal the candidate's 132 rows to the WRONG targets and re-run every endpoint,
standalone and blended.  If the permuted arm reproduces the gain, the gain is a gene-level
marginal and carries no per-target direction information.
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
    base = score(s, inc)
    B0 = np.load(os.path.join(HERE, "work", "coef_k0.npz"),
                 allow_pickle=True)["beta"].astype(np.float64)
    C = match(frame(dbl(B0)))                       # the +beta arm

    print(HDR)
    print(ROW.format(name="incumbent", **base))
    print(ROW.format(name="+beta k=0 standalone", **score(s, C)))

    n = len(s.perts)
    keys = ("sign_acc", "head_purity", "reach", "pds", "nmae")
    for f in (None, 0.05, 0.10):
        real = score(s, C if f is None else (1 - f) * inc + f * C)
        acc = {k: [] for k in keys}
        for seed in range(50):
            pm = np.random.default_rng(seed).permutation(n)
            P = C[pm]
            r = score(s, P if f is None else (1 - f) * inc + f * P)
            for k in keys:
                acc[k].append(r[k])
        lab = "standalone" if f is None else f"blend f={f:.2f}"
        print()
        print(f"{lab}: candidate vs 50 row-permutations of itself")
        print(f"  {'endpoint':<12} {'real':>9} {'permuted mean':>15} {'sd':>8} {'z':>7}")
        for k in keys:
            m, sd = float(np.mean(acc[k])), float(np.std(acc[k]))
            print(f"  {k:<12} {real[k]:9.4f} {m:15.4f} {sd:8.4f} "
                  f"{(real[k]-m)/max(sd,1e-9):+7.2f}")

    print()
    print("mechanism check: where does the +beta ordering put its head?")
    idx = np.argsort(-np.abs(C[0]))
    for p in (0,):
        pass
    heavy = np.zeros(0, dtype=int)
    downfrac = []
    for p in range(n):
        i = np.flatnonzero(s.real_sig[p] & s.adj[p])
        if i.size < 10:
            continue
        m = i[np.abs(C[p, i]) >= np.percentile(np.abs(C[p, i]), 90)]
        downfrac.append((s.real_lfc[p, m] < 0).mean())
        heavy = np.concatenate([heavy, m])
    allf = []
    for p in range(n):
        i = np.flatnonzero(s.real_sig[p] & s.adj[p])
        if i.size >= 10:
            allf.append((s.real_lfc[p, i] < 0).mean())
    print(f"  reference pool is {100*np.mean(allf):.1f}% down overall, but "
          f"{100*np.mean(downfrac):.1f}% down inside +beta's top decile")
    print(f"  +beta predicts down for "
          f"{100*(C[np.arange(n)[:, None], np.arange(C.shape[1])] < 0).mean():.1f}% of all "
          f"entries")


if __name__ == "__main__":
    main()
