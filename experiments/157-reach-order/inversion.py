"""157 Q1 — is the exp-094-era ORDER INVERSION still present under the shipped emitter?

094 measured, on the pre-a8 (minimum-variance / phi=0) emission:
    sign accuracy inside our call set  0.3446
    sign accuracy outside it           0.8259   (paired t = -54.0, inverted in 99.2%)
    Spearman(our p_adj, sign correct) +0.5009   (higher p_adj => MORE often correct)

Since then a8/a15 ship gamma-Poisson phi=4 at spread 0.25.  This recomputes the same three
statistics on any DE table, so a phi=0 table and a phi=4 table can be compared directly.

Deterministic given the table.  No emitter run here, no cell-eval2, no build, no submission.
"""
from __future__ import annotations
import os
import sys

import numpy as np
from scipy import stats

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
from debench import Surface  # noqa: E402

TABLE = os.environ.get("PRED_TABLE", f"{ROOT}/experiments/095-clean-harness/de_thinned.npz")
NS = [None, 68, 132, 240, 400]


def report(s: Surface, tag: str) -> None:
    P = len(s.perts)
    ins_n = ins_k = out_n = out_k = 0
    per_pert, rhos = [], []
    for p in range(P):
        m = s.real_sig[p] & s.adj[p]
        i = np.flatnonzero(m)
        if i.size == 0:
            continue
        ok = (np.sign(s.pred_lfc[p, i]) == np.sign(s.real_lfc[p, i]))
        sig = s.pred_sig[p, i]
        ins_n += int(sig.sum()); ins_k += int(ok[sig].sum())
        out_n += int((~sig).sum()); out_k += int(ok[~sig].sum())
        if sig.sum() >= 5 and (~sig).sum() >= 5:
            per_pert.append((ok[sig].mean(), ok[~sig].mean()))
        # Spearman(p_adj, correct) INSIDE the pool: positive => higher p_adj more often right
        q = s.pred_padj[p, i]
        good = np.isfinite(q)
        if good.sum() >= 20 and 0 < ok[good].sum() < good.sum():
            rhos.append(stats.spearmanr(q[good], ok[good].astype(float)).statistic)
    a = np.array(per_pert)
    t = stats.ttest_rel(a[:, 0], a[:, 1]) if a.size else None
    rho = np.array([r for r in rhos if np.isfinite(r)])
    print(f"\n=== {tag} ===")
    print(f"  pool (ref-significant, adjudicable) pairs      {ins_n + out_n:>8d}")
    print(f"  sign accuracy INSIDE our call set              {ins_k / max(ins_n,1):>8.4f}  (n={ins_n})")
    print(f"  sign accuracy OUTSIDE our call set             {out_k / max(out_n,1):>8.4f}  (n={out_n})")
    if t is not None:
        print(f"  paired t over {a.shape[0]:d} perturbations{'':>13}{t.statistic:>+8.1f}  p={t.pvalue:.2e}")
        print(f"  fraction of perturbations INVERTED (in<out)    {float((a[:,0] < a[:,1]).mean()):>8.4f}")
    print(f"  mean per-pert Spearman(pred p_adj, correct)    {rho.mean():>+8.4f}  "
          f"(median {np.median(rho):+.4f}, n={rho.size})")
    print(f"  fraction of perturbations with rho > 0         {float((rho > 0).mean()):>8.4f}")
    # pooled Spearman, matching 094's headline framing
    m = s.real_sig & s.adj
    q = np.where(np.isfinite(s.pred_padj), s.pred_padj, 1.0)[m]
    ok = (np.sign(s.pred_lfc) == np.sign(s.real_lfc))[m].astype(float)
    print(f"  POOLED Spearman(pred p_adj, correct)           "
          f"{stats.spearmanr(q, ok).statistic:>+8.4f}   [094 reported +0.5009]")


def main() -> None:
    print(f"table: {TABLE}")
    for N in NS:
        s = Surface(TABLE, n_real_target=N)
        nr, npd = s.counts()
        report(s, f"n_real = {'native BH' if N is None else N} "
                  f"(median n_real {np.median(nr):.0f}, median n_pred {np.median(npd):.0f})")


if __name__ == "__main__":
    main()
