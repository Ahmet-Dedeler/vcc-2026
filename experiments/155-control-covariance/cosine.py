"""155 — the powerful version of the question.

`plus.py` leaves one live thread: at a 5% blend the `+beta` arm moves `pds` +0.0055 while its
own row-permutation moves it -0.0102, i.e. a +2.2 sigma target-specific component.  `pds` is a
RANK statistic over 132 candidates and throws away most of the information; the underlying
quantity is simply whether a target's coefficient row points at that target's true effect.

So measure it directly: cosine(candidate row p, reference lfc row p) against the 131*132
off-diagonal cosines.  This is the same statistic AGENTS.md 4.2 calls the transfer cosine, it
is what `pds` is a coarse rank-summary of, and a paired test on it is far better powered than
50 permutations of a rank.

Every arm is scored on the reference's DE-significant genes only (where the reference has an
opinion) and on all scored genes, because the two can disagree.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bench import load_frame, score  # noqa: E402

s, inc, gdst, pok = load_frame()
n = len(s.perts)


def frame(M):
    P = np.zeros_like(s.pred_lfc)
    P[:, gdst] = M[:, gdst]
    return P


def dbl(M):
    M = M - M.mean(0, keepdims=True)
    return M - M.mean(1, keepdims=True)


def cos_table(P, cols):
    A = P[:, cols].astype(np.float64)
    B = s.real_lfc[:, cols].astype(np.float64)
    A /= np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-30)
    B /= np.maximum(np.linalg.norm(B, axis=1, keepdims=True), 1e-30)
    return A @ B.T


def report(name, P):
    # panel targets are removed from pds; do the same here
    sidx = set(int(i) for i in s.self_idx if i >= 0)
    cols = np.array([i for i in gdst if i not in sidx])
    S = cos_table(P, cols)
    d = np.diag(S)
    off = S[~np.eye(n, dtype=bool)]
    # paired: each target's own cosine minus the mean cosine of the other 131 predictions
    # against that same truth -- removes any per-target scale or magnitude effect
    paired = d - (S.sum(0) - d) / (n - 1)
    t = paired.mean() / (paired.std(ddof=1) / np.sqrt(n))
    print(f"  {name:<34} own {d.mean():+.5f}   off {off.mean():+.5f}   "
          f"paired {paired.mean():+.5f}   t({n-1}) = {t:+.2f}")
    return paired


def main() -> None:
    print("transfer cosine on the reference's full scored gene axis (panel targets removed)")
    B0 = np.load(os.path.join(HERE, "work", "coef_k0.npz"),
                 allow_pickle=True)["beta"].astype(np.float64)
    B20 = np.load(os.path.join(HERE, "work", "coef_k20.npz"),
                  allow_pickle=True)["beta"].astype(np.float64)
    report("INCUMBENT 104 baseline", inc)
    report("+beta k=0  (double-centered)", frame(dbl(B0)))
    report("-beta k=0  (the hypothesis)", frame(-dbl(B0)))
    report("+beta k=20 (double-centered)", frame(dbl(B20)))
    report("-beta k=20 (the hypothesis)", frame(-dbl(B20)))

    print()
    print("tightened permutation null for the 5% and 10% blends (500 seeds, +beta k=0)")
    a = np.linalg.norm(inc[:, gdst], axis=1)
    C = frame(dbl(B0))
    C *= (a / np.maximum(np.linalg.norm(C[:, gdst], axis=1), 1e-30))[:, None]
    for f in (0.05, 0.10):
        real = score(s, (1 - f) * inc + f * C)
        acc = {k: [] for k in ("reach", "pds", "head_purity", "sign_acc")}
        for seed in range(500):
            pm = np.random.default_rng(seed).permutation(n)
            r = score(s, (1 - f) * inc + f * C[pm])
            for k in acc:
                acc[k].append(r[k])
        print(f"  f = {f:.2f}")
        for k in ("sign_acc", "head_purity", "reach", "pds"):
            m, sd = float(np.mean(acc[k])), float(np.std(acc[k], ddof=1))
            print(f"    {k:<12} real {real[k]:.4f}   perm {m:.4f} +- {sd:.4f}   "
                  f"z = {(real[k]-m)/sd:+.2f}")



def stratified() -> None:
    """The estimator has no power on a 4-CPM target, so the a-priori-strongest test is the
    paired transfer cosine restricted to the best-measured targets."""
    cpm = np.load(os.path.join(HERE, "work", "coef_k20.npz"),
                  allow_pickle=True)["target_cpm"]
    order = np.argsort(-cpm)
    sidx = set(int(i) for i in s.self_idx if i >= 0)
    cols = np.array([i for i in gdst if i not in sidx])
    print()
    print("paired transfer cosine, stratified by how well the target itself is measured")
    print(f"{'arm':<22} {'stratum':<18} {'n':>4} {'paired cos':>12} {'se':>9} {'t':>7}")
    for lab, M in (("+beta k=0", None), ("+beta k=20", None)):
        k = 0 if "k=0" in lab else 20
        B = np.load(os.path.join(HERE, "work", f"coef_k{k}.npz"),
                    allow_pickle=True)["beta"].astype(np.float64)
        S = cos_table(frame(dbl(B)), cols)
        for lo, hi, nm in ((0, 33, "top quartile CPM"), (0, 66, "top half CPM"),
                           (0, 132, "all targets")):
            r = order[lo:hi]
            Ss = S[np.ix_(r, r)]
            d = np.diag(Ss)
            m = Ss.shape[0]
            paired = d - (Ss.sum(0) - d) / (m - 1)
            se = paired.std(ddof=1) / np.sqrt(m)
            print(f"{lab:<22} {nm:<18} {m:>4} {paired.mean():+12.5f} {se:9.5f} "
                  f"{paired.mean()/se:+7.2f}")
    S = cos_table(inc, cols)
    d = np.diag(S)
    paired = d - (S.sum(0) - d) / (n - 1)
    print(f"{'INCUMBENT':<22} {'all targets':<18} {n:>4} {paired.mean():+12.5f} "
          f"{paired.std(ddof=1)/np.sqrt(n):9.5f} "
          f"{paired.mean()/(paired.std(ddof=1)/np.sqrt(n)):+7.2f}")



if __name__ == "__main__":
    import sys as _s
    if "--strat" not in _s.argv:
        main()
    stratified()
