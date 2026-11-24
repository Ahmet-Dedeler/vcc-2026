"""165 — a higher-power reading of the same tables: reach on the SEED-AVERAGED order, and
the toll split into its two channels.

WHY.  `reso.py` measures that a single emission's `reach` on 85 perturbations carries a
clustered standard error of ~0.0037, which is the size of the whole effect: a23 itself (`mn`)
reads -0.00375 +/- 0.00373, t = -1.0.  Averaging `reach` over seeds shrinks that slowly,
because `k*` is a step function of one realised ordering.

Two things with far more power come out of the same data:

1. SEED-AVERAGED ORDER.  For each arm, average each gene's POSITION in the production
   confidence order (significant first, then ascending p_adj, ...) over all its seeds, and
   score `reach` through that averaged order.  The averaged order is a much better estimate of
   the arm's EXPECTED confidence ranking than any single draw, and `reach` computed on it is
   deterministic.  It is a SURROGATE -- `reach(E[order]) != E[reach(order)]` -- so it is
   gated on a known-answer row (a8 must still beat a15) and never converted to board units on
   its own.

2. CHANNEL SPLIT.  The mass correction can hurt `reach` two ways: by changing the ORDER genes
   are walked in, or by flipping the SIGN of a predicted lfc inside the pool.  Holding one of
   the two at a15's value isolates the other.  156 measured that the sign flips a23 causes are
   68.9% CORRECT, so the sign channel should HELP; if the split says otherwise, 156's claim is
   wrong.

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
from debench import Surface                                   # noqa: E402

ORDER = ["a15", "a8", "mn", "mn_half", "cst", "hp_p500", "hp_p2000", "hp_p6000",
         "hp_p9000", "hp_w2000", "eps_keep", "eps_tail", "inveff", "eff", "mn_c10"]


def positions(s):
    """Position of every gene in the production confidence order, per perturbation."""
    P, G = s.pred_p.shape
    out = np.empty((P, G), dtype=np.int32)
    for p in range(P):
        o = np.lexsort((s.genes, -np.abs(s.pred_lfc[p]), s.pred_p[p],
                        s.pred_padj[p], ~s.pred_sig[p]))
        out[p, o] = np.arange(G, dtype=np.int32)
    return out


def per_pert_reach(s, key=None, pred_lfc=None):
    lfc = s.pred_lfc if pred_lfc is None else pred_lfc
    out = np.full(len(s.perts), np.nan)
    for p in range(len(s.perts)):
        pool = np.flatnonzero(s.real_sig[p])
        if pool.size == 0:
            continue
        pool = pool[s._order(p, pool, key, None)]
        pool = pool[s.adj[p, pool]]
        if pool.size == 0:
            out[p] = 0.0
            continue
        ok = np.sign(lfc[p, pool]) == np.sign(s.real_lfc[p, pool])
        pur = np.cumsum(ok) / np.arange(1, pool.size + 1)
        deep = np.flatnonzero(pur >= 0.9)
        out[p] = (int(deep[-1] + 1) if deep.size else 0) / int(s.real_sig[p].sum())
    return out


def main() -> None:
    seeds = {a: sorted(int(re.search(r"\.s(\d+)\.npz$", f).group(1))
                       for f in glob.glob(f"{WORK}/de.{a}.s*.npz"))
             for a in ORDER if glob.glob(f"{WORK}/de.{a}.s*.npz")}
    arms = [a for a in ORDER if a in seeds]
    print("seeds per arm: " + ", ".join(f"{a}:{len(seeds[a])}" for a in arms))

    S0 = {a: Surface(f"{WORK}/de.{a}.s{seeds[a][0]}.npz") for a in arms}

    def averaged(arm, seedlist):
        acc = lacc = None
        for s in seedlist:
            S = Surface(f"{WORK}/de.{arm}.s{s}.npz")
            pz = positions(S).astype(np.float64)
            acc = pz if acc is None else acc + pz
            lacc = S.pred_lfc if lacc is None else lacc + S.pred_lfc
        return acc / len(seedlist), lacc / len(seedlist)

    mpos, mlfc = {}, {}
    for a in arms:
        mpos[a], mlfc[a] = averaged(a, seeds[a])
    # THE SURROGATE'S LEVEL DEPENDS ON HOW MANY SEEDS WENT INTO THE AVERAGE: more seeds ->
    # a better-estimated order -> higher reach (a15 reads 0.17566 at 2 seeds and 0.18303 at
    # 12). So every arm is compared against an a15 baseline built from the SAME seed list.
    # Comparing a 2-seed arm against a 12-seed a15 manufactures a large negative delta out of
    # nothing, which is what an earlier pass of this script did.
    a15_by_n = {}
    for a in arms:
        k = tuple(seeds[a])
        if k not in a15_by_n:
            a15_by_n[k] = averaged("a15", seeds[a])

    base = S0["a15"]
    print("\n" + "=" * 100)
    print("(A) reach on the SEED-AVERAGED confidence order, arm's own averaged pred_lfc")
    print("=" * 100)
    print(f"{'arm':<11}{'seeds':>7}{'reach(mean order)':>19}{'a15 same n':>13}"
          f"{'d vs a15':>11}{'clustered se':>14}{'t':>7}")
    res = {}
    R0 = {k: per_pert_reach(base, key=-v[0], pred_lfc=v[1]) for k, v in a15_by_n.items()}
    for a in arms:
        S = S0[a]
        S.real_sig = base.real_sig          # identical by construction; make it explicit
        S.real_lfc = base.real_lfc
        S.adj = base.adj
        r0 = R0[tuple(seeds[a])]
        r = per_pert_reach(S, key=-mpos[a], pred_lfc=mlfc[a])
        d = r - r0
        fin = np.isfinite(d)
        dd = float(np.mean(d[fin]))
        se = float(np.std(d[fin], ddof=1) / np.sqrt(fin.sum()))
        res[a] = (float(np.nanmean(r)), dd, se)
        print(f"{a:<11}{len(seeds[a]):>7}{np.nanmean(r):>19.5f}"
              f"{np.nanmean(r0):>13.5f}"
              + (f"{'—':>11}{'—':>14}{'—':>7}" if a == "a15" else
                 f"{dd:>+11.5f}{se:>14.5f}{dd / se if se else 0:>7.2f}"))
    if "a8" in res:
        ok = res["a8"][1] > 0
        print(f"\n  KNOWN-ANSWER: a8 must still beat a15 on this surrogate — "
              f"{res['a8'][1]:+.5f}  {'PASS' if ok else 'FAIL'}")
        if not ok:
            print("  -> the surrogate does not carry the one contrast we know the sign of; "
                  "do not use it.")

    print("\n" + "=" * 100)
    print("(B) channel split, on the seed-averaged order: ORDER only vs SIGN only")
    print("=" * 100)
    print(f"{'arm':<11}{'both':>11}{'order only':>12}{'sign only':>11}"
          f"{'order+sign':>12}")
    for a in arms:
        if a == "a15":
            continue
        S = S0[a]
        S.real_sig, S.real_lfc, S.adj = base.real_sig, base.real_lfc, base.adj
        mp0, ml0 = a15_by_n[tuple(seeds[a])]
        both = np.nanmean(per_pert_reach(S, key=-mpos[a], pred_lfc=mlfc[a]))
        ordo = np.nanmean(per_pert_reach(S, key=-mpos[a], pred_lfc=ml0))
        sgn = np.nanmean(per_pert_reach(S, key=-mp0, pred_lfc=mlfc[a]))
        b0 = np.nanmean(R0[tuple(seeds[a])])
        print(f"{a:<11}{both - b0:>+11.5f}{ordo - b0:>+12.5f}{sgn - b0:>+11.5f}"
              f"{(ordo - b0) + (sgn - b0):>+12.5f}")
    print("\n  'order only' holds pred_lfc at a15's and moves only the walk order;")
    print("  'sign only' holds the walk order at a15's and moves only the predicted signs.")


if __name__ == "__main__":
    main()
