"""155 — close out the one cell that looked positive.

`blend.py` §2 reported the k=0 candidate at sign accuracy 0.5487 on the 33 best-measured
targets against a permutation null of 0.4507 (z = +2.88).  That null is DEFLATED: it deals
rows from all 132 targets into a stratum of 33, and a target's coefficient row carries a
global offset that depends on how well the target is measured, so the permuted rows are
systematically worse for reasons that have nothing to do with target identity.

The matched null permutes WITHIN the stratum, where every row has a comparable offset.  Also
reported: the sign composition of the candidate, which decides whether anything here is more
than the reference pool's own 60% down-regulation marginal.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bench import head_purity, load_frame  # noqa: E402

s, inc, gdst, pok = load_frame()
RS, AD = s.real_sig.copy(), s.adj.copy()


def frame(M):
    P = np.zeros_like(s.pred_lfc)
    P[:, gdst] = M[:, gdst]
    return P


def dbl(M):
    M = M - M.mean(0, keepdims=True)
    return M - M.mean(1, keepdims=True)


def restrict(rows):
    keep = np.zeros_like(RS)
    keep[rows] = True
    s.real_sig, s.adj = RS & keep, AD & keep


def q(P):
    s.pred_lfc = P
    return s.sign_accuracy(), head_purity(s, P)[0]


def main() -> None:
    cpm = np.load(os.path.join(HERE, "work", "coef_k20.npz"),
                  allow_pickle=True)["target_cpm"]
    order = np.argsort(-cpm)
    arms = {}
    for k in (0, 20):
        B = np.load(os.path.join(HERE, "work", f"coef_k{k}.npz"),
                    allow_pickle=True)["beta"].astype(np.float64)
        arms[k] = frame(-dbl(B))

    print("WITHIN-STRATUM permutation null (rows dealt only among the stratum's own "
          "targets, 200 seeds)")
    print(f"{'stratum':<22} {'k':>3} {'cand acc':>9} {'null acc':>16} {'z':>6} "
          f"{'cand head':>10} {'null head':>16} {'z':>6} {'marginal':>9}")
    for lo, hi, name in ((0, 33, "top quartile CPM"), (0, 66, "top half CPM"),
                         (0, 132, "all 132 targets")):
        rows = order[lo:hi]
        restrict(rows)
        pool = np.concatenate([s.real_lfc[p][s.real_sig[p] & s.adj[p]] for p in rows])
        marg = max((pool < 0).mean(), (pool > 0).mean())
        for k in (0, 20):
            P = arms[k]
            a, h = q(P)
            na, nh = [], []
            for seed in range(200):
                pm = np.arange(len(s.perts))
                pm[rows] = rows[np.random.default_rng(seed).permutation(rows.size)]
                aa, hh = q(P[pm])
                na.append(aa)
                nh.append(hh)
            za = (a - np.mean(na)) / max(np.std(na), 1e-9)
            zh = (h - np.mean(nh)) / max(np.std(nh), 1e-9)
            print(f"{name:<22} {k:>3} {a:9.4f} {np.mean(na):9.4f}+-{np.std(na):.4f} "
                  f"{za:+6.2f} {h:10.4f} {np.mean(nh):9.4f}+-{np.std(nh):.4f} {zh:+6.2f} "
                  f"{marg:9.4f}")
    s.real_sig, s.adj = RS, AD

    print()
    print("sign composition on the reference DE pool -- is the candidate anything but the "
          "marginal?")
    for k in (0, 20):
        P = arms[k]
        neg = num = 0
        for p in range(len(s.perts)):
            i = np.flatnonzero(RS[p] & AD[p])
            neg += int((P[p, i] < 0).sum())
            num += i.size
        print(f"  k={k:<3} candidate predicts {100*neg/num:5.1f}% down; "
              f"reference pool is 59.7% down")
    negi = num = 0
    for p in range(len(s.perts)):
        i = np.flatnonzero(RS[p] & AD[p])
        negi += int((inc[p, i] < 0).sum())
        num += i.size
    print(f"  incumbent predicts {100*negi/num:5.1f}% down")


if __name__ == "__main__":
    main()
