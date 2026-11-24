"""165 — can this instrument resolve a PLACEMENT difference in `reach` at all?

`price.py` reports arm-to-arm raw `reach` deltas of 0.002-0.005 against a seed-to-seed spread
of 0.009.  That is the wrong comparison (the seed spread is a common shift both arms share),
but it is the right warning: `reach = mean_p(k*_p / n_real_p)` and `k*` is the DEEPEST prefix
at >=90% purity, a step function of the ordering.  One gene moving across the boundary can
change a perturbation's contribution by a large fraction, and 157 §3a measured that 98% of our
reach comes from the small-`n_real` tercile — so the effective sample size is far below 85.

This script does three things, none of which need another emission:

1. PAIRED, CLUSTERED ERROR BARS.  The arm delta is a mean over perturbations of a paired
   quantity `k*_p(arm)/n_real_p - k*_p(a15)/n_real_p` at the SAME seed.  The standard error is
   computed by clustering on perturbation (each perturbation's per-seed deltas averaged first),
   which is the honest error bar for a mean over a heterogeneous, seed-repeated panel.

2. PROPORTIONAL POWER MATCHING (157 §3b's instrument, `PowerSurface`).  `f` scales each
   perturbation's reference pool by a common factor, preserving the heterogeneity `reach`
   averages over.  `f = 1` is the native BH set — the known-answer row.  `f = 3.5` is the
   challenge-like power regime, where the walk is longer and `k*` is less of a step function.

3. A SMOOTH SURROGATE, `ordAUC`.  Within each perturbation's reference DE pool, the probability
   that a sign-CORRECT gene is ranked ahead of a sign-INCORRECT one in our confidence order.
   It is the ordering content `reach` is trying to reward, without the >=90%-prefix step, so it
   has far more power.  It is NOT the scored metric and is never converted to board units — it
   is used only to decide whether a placement difference exists at all.

Reads only the tables `arms.py` wrote.  No emitter, no cell-eval2, no build, no submission.
"""
from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
sys.path.insert(0, f"{ROOT}/experiments/157-reach-order")
from debench import Surface                                   # noqa: E402
from power import PowerSurface                                # noqa: E402

ORDER = ["a15", "a8", "mn", "mn_half", "cst", "hp_p500", "hp_p2000", "hp_p6000",
         "hp_p9000", "hp_w2000", "eps_keep", "eps_tail", "inveff", "eff", "mn_c10"]
FACTORS = (1.0, 3.5)


def per_pert_reach(s):
    """`k*_p / n_real_p` for every perturbation, NaN where the pool is empty."""
    out = np.full(len(s.perts), np.nan)
    for p in range(len(s.perts)):
        pool = np.flatnonzero(s.real_sig[p])
        if pool.size == 0:
            continue
        pool = pool[s._order(p, pool, None, None)]
        pool = pool[s.adj[p, pool]]
        if pool.size == 0:
            out[p] = 0.0
            continue
        ok = np.sign(s.pred_lfc[p, pool]) == np.sign(s.real_lfc[p, pool])
        pur = np.cumsum(ok) / np.arange(1, pool.size + 1)
        deep = np.flatnonzero(pur >= 0.9)
        out[p] = (int(deep[-1] + 1) if deep.size else 0) / int(s.real_sig[p].sum())
    return out


def per_pert_ordauc(s):
    """P(a sign-correct pool gene ranks ahead of a sign-incorrect one), per perturbation."""
    out = np.full(len(s.perts), np.nan)
    for p in range(len(s.perts)):
        pool = np.flatnonzero(s.real_sig[p])
        if pool.size == 0:
            continue
        pool = pool[s._order(p, pool, None, None)]
        pool = pool[s.adj[p, pool]]
        ok = np.sign(s.pred_lfc[p, pool]) == np.sign(s.real_lfc[p, pool])
        n1, n0 = int(ok.sum()), int((~ok).sum())
        if n1 == 0 or n0 == 0:
            continue
        # rank 0 = most confident; AUC over the ordering, ties impossible (total order)
        r = np.arange(pool.size)
        out[p] = 1.0 - (r[ok].sum() - n1 * (n1 - 1) / 2) / (n1 * n0)
    return out


def report(name, vals, arms, base="a15"):
    """`vals[(arm, seed)]` -> per-perturbation vector. Paired, clustered on perturbation."""
    seeds = sorted({s for (_, s) in vals})
    print(f"\n{'arm':<11}{'level':>10}{'d vs a15':>11}{'clustered se':>14}"
          f"{'t':>7}{'per-seed sd':>13}{'n seeds':>9}")
    res = {}
    for a in arms:
        lv = np.nanmean([np.nanmean(vals[(a, s)]) for s in seeds])
        if a == base:
            print(f"{a:<11}{lv:>10.5f}{'—':>11}{'—':>14}{'—':>7}"
                  f"{np.std([np.nanmean(vals[(a, s)]) for s in seeds]):>13.5f}"
                  f"{len(seeds):>9}")
            continue
        D = np.stack([vals[(a, s)] - vals[(base, s)] for s in seeds])   # (S, P)
        per_pert = np.nanmean(D, axis=0)
        fin = np.isfinite(per_pert)
        d = float(np.nanmean(per_pert[fin]))
        se = float(np.nanstd(per_pert[fin], ddof=1) / np.sqrt(fin.sum()))
        psd = float(np.std([np.nanmean(D[i]) for i in range(D.shape[0])], ddof=1)) \
            if len(seeds) > 1 else float("nan")
        res[a] = (d, se)
        print(f"{a:<11}{lv:>10.5f}{d:>+11.5f}{se:>14.5f}{d / se if se else 0:>7.2f}"
              f"{psd:>13.5f}{len(seeds):>9}")
    return res


def main() -> None:
    seeds = sorted({int(re.search(r"\.s(\d+)\.npz$", f).group(1))
                    for f in glob.glob(f"{WORK}/de.a15.s*.npz")})
    arms = [a for a in ORDER
            if all(os.path.exists(f"{WORK}/de.{a}.s{s}.npz") for s in seeds)]
    partial = [a for a in ORDER if a not in arms
               and glob.glob(f"{WORK}/de.{a}.s*.npz")]
    print(f"seeds present for a15: {seeds}")
    print(f"arms present at ALL those seeds: {arms}")
    if partial:
        sub = sorted({int(re.search(r"\.s(\d+)\.npz$", f).group(1))
                      for f in glob.glob(f"{WORK}/de.{partial[0]}.s*.npz")})
        print(f"arms present only at {sub}: {partial}  -> reported separately below")

    for f in FACTORS:
        print("\n" + "=" * 104)
        print(f"REACH, proportional power matching f = {f}"
              + ("   (native BH set — the known-answer row)" if f == 1.0
                 else "   (challenge-like power)"))
        print("=" * 104)
        vals = {}
        for a in arms:
            for s in seeds:
                S = PowerSurface(f"{WORK}/de.{a}.s{s}.npz", factor=f) if f != 1.0 \
                    else Surface(f"{WORK}/de.{a}.s{s}.npz")
                vals[(a, s)] = per_pert_reach(S)
        report(f"reach f={f}", vals, arms)
        if f == 1.0 and "mn" in arms:
            # The decision is head-protected vs MINNORM, not vs a15: a15 is the arm we
            # already know we are leaving. Same paired, clustered construction.
            print("\n  ... the same table with `mn` (= a23) as the base:")
            report(f"reach f={f} vs mn", vals, [a for a in arms if a != "a15"], base="mn")

    print("\n" + "=" * 104)
    print("ordAUC — the smooth ordering surrogate (NOT a scored metric, never converted)")
    print("=" * 104)
    vals = {}
    for a in arms:
        for s in seeds:
            vals[(a, s)] = per_pert_ordauc(Surface(f"{WORK}/de.{a}.s{s}.npz"))
    report("ordAUC", vals, arms)

    if partial:
        sub = sorted({int(re.search(r"\.s(\d+)\.npz$", f).group(1))
                      for f in glob.glob(f"{WORK}/de.{partial[0]}.s*.npz")})
        print("\n" + "=" * 104)
        print(f"the same two tables restricted to the {len(sub)}-seed subset "
              f"({partial + ['a15']})")
        print("=" * 104)
        arms2 = ["a15"] + partial
        for f in FACTORS:
            vals = {}
            for a in arms2:
                for s in sub:
                    S = PowerSurface(f"{WORK}/de.{a}.s{s}.npz", factor=f) if f != 1.0 \
                        else Surface(f"{WORK}/de.{a}.s{s}.npz")
                    vals[(a, s)] = per_pert_reach(S)
            print(f"\nREACH f = {f}")
            report(f"reach f={f}", vals, arms2)
        vals = {}
        for a in arms2:
            for s in sub:
                vals[(a, s)] = per_pert_ordauc(Surface(f"{WORK}/de.{a}.s{s}.npz"))
        print("\nordAUC")
        report("ordAUC", vals, arms2)


if __name__ == "__main__":
    main()
