"""165 — what actually populates the `reach` head, and where the mass correction lands.

`A23-BOARD-RESULT.md` states the mechanism as: minimum-L2 puts `eps ∝ -exp(b_ctrl)` on the
HIGHEST-EXPRESSED genes, those carry the most counts, so they have the smallest emitted
p-values, so they populate the head `reach` reads.  157 §5a states a different mechanism for
the same head: under a15's shipped tier policy the ~1,900-gene `spread 1.0 / phi 4` TAIL block
is significant at zero predicted effect with median |z| = 7.4 against the keep block's 1.7, and
supplies 79.7% of the top ten.

Those two stories predict different head-protection rules, so this script measures which is
true on the surface that reproduces the board's `reach`, and writes `work/head.npz` (the
per-target confidence order `reach` walks) for `price.py`.

Reads only the tables `arms.py` wrote.  No emitter, no cell-eval2, no build, no submission.
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
from debench import Surface                                   # noqa: E402

TS = 5.0e4
SEED = 90_000


def main() -> None:
    s = Surface(f"{WORK}/de.a15.s{SEED}.npz")
    tab = np.load(f"{WORK}/de.a15.s{SEED}.npz", allow_pickle=True)
    spread = tab["spread"]                       # (T, G) on the scored axis
    z = np.load(f"{WORK}/profiles_all.npz", allow_pickle=True)
    b = z["b_ctrl"]
    keep = z["keep"]
    kidx = np.flatnonzero(keep)
    p_a15 = z["p_a15"].astype(np.float64)
    gate = z["gate"]
    T, G = s.pred_p.shape
    print(f"h1-test: {T} perturbations x {G} scored genes; "
          f"n_real median {np.median([int(s.real_sig[p].sum()) for p in range(T)]):.0f}")

    # ---- 1. the head reach actually reads ---------------------------------------------
    print("\n" + "=" * 96)
    print("(1) what populates the head of `reach`'s walk (a15, the shipped emission)")
    print("=" * 96)
    w_rank = np.empty((T, G))                    # percentile of control expression
    expr = b[keep]
    er = np.argsort(np.argsort(expr)) / (G - 1)  # 1.0 = highest expressed
    rows = {k: [] for k in ("tail", "expr", "kstar", "n_real")}
    depth_expr, depth_tail = [], []
    for p in range(T):
        pool = np.flatnonzero(s.real_sig[p])
        if pool.size == 0:
            continue
        order = s._order(p, pool, None, None)
        pool = pool[order]
        pool = pool[s.adj[p, pool]]
        if pool.size == 0:
            continue
        ok = np.sign(s.pred_lfc[p, pool]) == np.sign(s.real_lfc[p, pool])
        pur = np.cumsum(ok) / np.arange(1, pool.size + 1)
        deep = np.flatnonzero(pur >= 0.9)
        k = int(deep[-1] + 1) if deep.size else 0
        head = pool[:max(k, 10)]
        rows["tail"].append(float((spread[p, head] > 0.5).mean()))
        rows["expr"].append(float(er[head].mean()))
        rows["kstar"].append(k)
        rows["n_real"].append(int(s.real_sig[p].sum()))
        depth_expr.append(er[pool[:10]].mean())
        depth_tail.append(float((spread[p, pool[:10]] > 0.5).mean()))
    print(f"  k* median {np.median(rows['kstar']):.0f} of n_real median "
          f"{np.median(rows['n_real']):.0f}   -> reach {np.mean(rows['kstar']):.2f}/"
          f"{np.mean(rows['n_real']):.2f}")
    print(f"  share of the pure head drawn from the spread-1.0 TAIL block: "
          f"{100 * np.mean(rows['tail']):.1f}%   (157 §5a measured 79.7% for the top ten on "
          f"repaired HCT116)")
    print(f"  share of the first TEN pool genes from the tail block: "
          f"{100 * np.mean(depth_tail):.1f}%")
    print(f"  mean control-expression percentile of the head: "
          f"{100 * np.mean(rows['expr']):.1f}%  (50% = no expression preference; "
          f"a23's stated mechanism needs this near 100%)")
    print(f"  ... of the first ten: {100 * np.mean(depth_expr):.1f}%")
    tailfrac = float((spread > 0.5).mean())
    print(f"  for reference the tail block is {100 * tailfrac:.1f}% of scored genes, so "
          f"'no preference' would read {100 * tailfrac:.1f}%")

    # ---- 2. is the global p-order driven by expression or by the spread block? ---------
    print("\n" + "=" * 96)
    print("(2) the GLOBAL emitted confidence order: what predicts it")
    print("=" * 96)
    from scipy.stats import spearmanr
    rho_expr, rho_spread, top_tail, top_expr = [], [], [], []
    for p in range(T):
        pv = s.pred_p[p]
        fin = np.isfinite(pv)
        rho_expr.append(spearmanr(pv[fin], expr[fin]).statistic)
        rho_spread.append(spearmanr(pv[fin], spread[p][fin]).statistic)
        o = np.lexsort((np.arange(G), pv))
        top_tail.append(float((spread[p, o[:2000]] > 0.5).mean()))
        top_expr.append(float(er[o[:2000]].mean()))
    print(f"  Spearman(pred p, control expression) mean {np.mean(rho_expr):+.4f}")
    print(f"  Spearman(pred p, emission spread)     mean {np.mean(rho_spread):+.4f}")
    print(f"  of the top 2000 by p: {100 * np.mean(top_tail):.1f}% are tail-block genes "
          f"(the tail block is only {100 * tailfrac:.1f}% of the axis)")
    print(f"  their mean expression percentile: {100 * np.mean(top_expr):.1f}%")

    # ---- 3. where minnorm's eps lands -------------------------------------------------
    print("\n" + "=" * 96)
    print("(3) where the correction lands, per placement, against that order")
    print("=" * 96)
    names = [k[2:] for k in z.files if k.startswith("p_")]
    order_full = np.zeros((T, G), dtype=np.int32)
    for p in range(T):
        order_full[p] = np.lexsort((np.arange(G), s.pred_p[p]))
    np.savez_compressed(f"{WORK}/head.npz", head=kidx[order_full], keep=keep)
    print(f"  wrote {WORK}/head.npz (per-target confidence order, full-axis indices)")
    print(f"\n{'arm':<11}{'|dp|/|p| %':>12}{'top-500 by p':>14}{'top-2000':>10}"
          f"{'tail block':>12}{'top-500 by expr':>17}")
    for a in names:
        if a in ("a15", "a8"):
            continue
        dp = np.abs(z[f"p_{a}"].astype(np.float64) - p_a15)[:, keep]
        rel = np.median(np.linalg.norm(dp, axis=1)
                        / np.linalg.norm(p_a15[:, keep], axis=1))
        s500 = s2000 = stail = sexpr = 0.0
        eo = np.argsort(-expr)
        for p in range(T):
            d = dp[p]
            tot = max(d.sum(), 1e-300)
            s500 += d[order_full[p][:500]].sum() / tot / T
            s2000 += d[order_full[p][:2000]].sum() / tot / T
            stail += d[spread[p] > 0.5].sum() / tot / T
            sexpr += d[eo[:500]].sum() / tot / T
        print(f"{a:<11}{100 * rel:>12.4f}{100 * s500:>13.1f}%{100 * s2000:>9.1f}%"
              f"{100 * stail:>11.1f}%{100 * sexpr:>16.1f}%")
    print(f"\n  'tail block' column baseline if the correction were spread uniformly over "
          f"the gate: {100 * tailfrac:.1f}%")


if __name__ == "__main__":
    main()
