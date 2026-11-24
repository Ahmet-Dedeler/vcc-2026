"""188 — the surface question, settled: `nmae`'s DE POOL COMPOSITION is the whole disagreement.

`identity188.py` produced two numbers that cannot both describe the same arm:

    h1-test, pool = reference BH<0.05          dU/D = **+0.0383**   d nmae = +0.04633
    PLOLO,   pool = top-2000 |lfc| per row     dU/D = **-0.0147**   d nmae = +0.00290

Opposite SIGNS on the L1 budget term.  Sharpening in the builder's `h` space moves magnitude
toward large-|h| genes; whether that ADDS L1 to the scored pool or REMOVES it depends entirely
on how the pool was selected.  A Wilcoxon BH set is selected by SIGNIFICANCE, which rises with
expression, and `|h| = log1p(TS*p_pert) - log1p(TS*p_ctrl)` also rises with expression at fixed
fold change.  A top-|lfc| set is selected the other way.

So one of the two pools is the challenge's and the other is mine.  This script decides it, in
two steps that do not depend on each other:

  A. On h1-test, where the reference's OWN BH<0.05 set is on disk, measure which construction it
     resembles: top-n by |lfc|, by |h|, or by a Wilcoxon-power proxy |lfc|*sqrt(cpm), at that
     row's own n.  This is a property of the scorer's DE call, not of any prediction.
  B. On PLOLO, re-price the arms under each pool construction at two sizes (89 = h1-test's
     median n_real, and 2,000 = 098 §4's estimate of the challenge's mean).  If the pool that
     wins step A reproduces h1-test's SIGN, the two surfaces never disagreed about the arm --
     only about which genes are scored -- and h1-test's reading is the one that transfers.

This is the question 169 §10 lists as "the open question, and it is the whole result".
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
sys.path.insert(0, HERE)
sys.path.insert(0, E169)
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
from debench import Surface                                    # noqa: E402
from budget188 import FOLDS, PDS_SPAN, deliver, load, pds_of    # noqa: E402
from mapfit import to_lfc                                       # noqa: E402
from sharp import sharpen, sharpen_lfc                           # noqa: E402

TS = 50_000.0


def lfc_to_h(L, cpm):
    base = TS * 1e-6 * np.asarray(cpm, dtype=np.float64)
    return np.log1p(base[None, :] * np.exp2(np.nan_to_num(L))) - np.log1p(base)[None, :]


def topn(key, n_per_row):
    out = np.zeros(key.shape, bool)
    for i in range(key.shape[0]):
        n = int(n_per_row[i])
        if n <= 0:
            continue
        n = min(n, key.shape[1])
        out[i, np.argpartition(-key[i], n - 1)[:n]] = True
    return out


def umd(P, R, pool):
    U, M, D = [], [], []
    for i in range(P.shape[0]):
        s = np.flatnonzero(pool[i])
        if s.size < 10:
            continue
        h, r = np.nan_to_num(P[i, s]), R[i, s]
        d = np.abs(r).sum()
        if not np.isfinite(d) or d == 0:
            continue
        m = np.where(np.sign(h) == np.sign(r), np.minimum(np.abs(h), np.abs(r)), 0.0)
        U.append(np.abs(h).sum()); M.append(m.sum()); D.append(d)
    return np.array(U), np.array(M), np.array(D)


def main() -> None:
    # ------------------------------------------------------------------ A. h1-test
    print("=" * 104)
    print("A.  WHAT KIND OF SET IS THE REFERENCE'S OWN BH<0.05 POOL?  (h1-test, scorer's "
          "call, no prediction)")
    E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
    pos = {str(g): i for i, g in enumerate(E["genes"])}
    s = Surface(os.path.join(E169, "de169.q1.00_mf.npz"))
    cpm = np.array([E["cpm"][pos[g]] for g in s.genes], dtype=np.float64)
    R = np.nan_to_num(s.real_lfc)
    Rh = lfc_to_h(R, cpm)
    sig = s.real_sig & np.isfinite(s.real_lfc)
    n = sig.sum(1)
    keys = dict(abs_lfc=np.abs(R), abs_h=np.abs(Rh),
                wald=np.abs(R) * np.sqrt(cpm)[None, :], cpm_only=np.broadcast_to(
                    cpm[None, :], R.shape).copy())
    print(f"    {len(s.perts)} perturbations; n_real per row: median {np.median(n):.0f}, "
          f"mean {n.mean():.0f}, min {n.min()}, max {n.max()}")
    print(f"    {'selection key':<14}{'overlap with the BH set':>26}")
    ovs = {}
    for k, key in keys.items():
        t = topn(key, n)
        o = float(np.mean([(t[i] & sig[i]).sum() / max(n[i], 1) for i in range(len(n))
                           if n[i] > 0]))
        ovs[k] = o
        print(f"    {k:<14}{o:>26.4f}")
    mc = float(np.median(cpm))
    inpool = float(np.median(np.concatenate([cpm[sig[i]] for i in range(len(n))
                                             if n[i] > 0])))
    print(f"    median control CPM: all genes {mc:.2f}, genes IN the BH pool {inpool:.2f}  "
          f"-> the pool is {inpool/mc:.1f}x more expressed than the gate")

    # ------------------------------------------------------------------ B. PLOLO
    print(f"\n{'=' * 104}")
    print("B.  PLOLO, the SAME arms, re-priced under each pool construction")
    z = np.load(os.path.join(WORK, "chalframes.npz"))
    cg = z["gA"] & z["gB"] & z["gC"]
    bA, gA = z["bA"], z["gA"]
    ARMS = ("q1.00", "q1.15", "q1.15L", "q1.25", "q1.10L", "q1.20L")
    POOLS = [("lfc", 89), ("lfc", 2000), ("h", 89), ("h", 2000),
             ("wald", 89), ("wald", 2000)]
    acc = {}
    for held in FOLDS:
        F = load(held)
        mk = np.zeros(cg.size, bool)
        mk[F["scored"]] = True
        c = np.flatnonzero(mk & cg)
        bh = F["b"]
        cpmf = 1e6 * np.expm1(bh[c]) / float(np.expm1(bh).sum())
        LR = to_lfc(F["truth"], bh, gA)[:, c]
        HR = F["truth"][:, c]
        pools = {}
        for kind, N in POOLS:
            key = dict(lfc=np.abs(LR), h=np.abs(HR),
                       wald=np.abs(LR) * np.sqrt(cpmf)[None, :])[kind]
            pools[(kind, N)] = topn(key, np.full(LR.shape[0], N))
        P0 = F["pred"]
        for name in ARMS:
            Q = {"q1.00": lambda: np.where(gA[None, :], P0, 0.0),
                 "q1.15": lambda: sharpen(P0, 1.15, gA, mode="row"),
                 "q1.25": lambda: sharpen(P0, 1.25, gA, mode="row"),
                 "q1.15L": lambda: sharpen_lfc(P0, 1.15, bA, gA, mode="row"),
                 "q1.10L": lambda: sharpen_lfc(P0, 1.10, bA, gA, mode="row"),
                 "q1.20L": lambda: sharpen_lfc(P0, 1.20, bA, gA, mode="row")}[name]()
            LP = to_lfc(deliver(Q, bA, gA), bA, gA)[:, c]
            for pk in pools:
                U, M, D = umd(LP, LR, pools[pk])
                acc.setdefault((name, pk), []).append((U, M, D))
        # pool composition diagnostic
        for pk in pools:
            m = pools[pk]
            print(f"    [{held}] pool {pk[0]}{pk[1]:<5} median ctrl CPM in pool "
                  f"{np.median(np.concatenate([cpmf[m[i]] for i in range(m.shape[0])])):8.2f} "
                  f"(gate median {np.median(cpmf):.2f})", flush=True)

    print(f"\n    {'pool':<12}{'arm':<9}{'nmae*':>9}{'d nmae*':>10}{'dU/D':>10}{'dM/D':>10}"
          f"{'e':>9}{'d nmae overall (span .60)':>27}")
    out = {}
    for kind, N in POOLS:
        b = acc[("q1.00", (kind, N))]
        U0 = np.concatenate([x[0] for x in b]); M0 = np.concatenate([x[1] for x in b])
        D0 = np.concatenate([x[2] for x in b])
        nm0 = float(np.mean(1 + (U0 - 2 * M0) / D0))
        for name in ARMS:
            v = acc[(name, (kind, N))]
            U = np.concatenate([x[0] for x in v]); M = np.concatenate([x[1] for x in v])
            nm = float(np.mean(1 + (U - 2 * M) / D0))
            du = float(np.mean((U - U0) / D0)); dm = float(np.mean((M - M0) / D0))
            out[f"{kind}{N}|{name}"] = dict(nmae=nm, dnmae=nm - nm0, dU=du, dM=dm)
            if name == "q1.00":
                continue
            print(f"    {kind + str(N):<12}{name:<9}{nm:>9.5f}{nm-nm0:>+10.5f}{du:>+10.5f}"
                  f"{dm:>+10.5f}{(dm/du if abs(du) > 1e-12 else float('nan')):>9.3f}"
                  f"{-(nm-nm0)/0.60/6:>+27.6f}")

    print(f"\n    h1-test's own reading, for comparison: q1.15 d nmae +0.04633 "
          f"(dU +0.03828, dM -0.00402), q1.15L +0.00855 (dU +0.00221, dM -0.00317)")
    with open(os.path.join(WORK, "pool188.json"), "w") as fh:
        json.dump(dict(h1_overlap=ovs, plolo=out), fh, indent=1, default=float)
    print(f"\nwrote {WORK}/pool188.json")


if __name__ == "__main__":
    sys.exit(main())
