"""172 — anatomy of `reach` on h1-test: where does the number physically come from?

Not an arm.  This exists because the brief's headline gap (0.235 -> 0.646) is only actionable
if you know WHICH targets and WHICH prefix lengths carry it, and because "precision at the top
is the right objective" needs a number attached to "the top".

`R_p = k*_p / n_real_p` with `k*` the deepest prefix at >=90% purity.  At a front accuracy near
0.62 a pure prefix is a handful of genes, so the mean is dominated by targets with SMALL pools
that happen to get their first few genes right.  This prints that distribution.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from harness172 import H1Test, check      # noqa: E402


def kstar(T, key):
    """Per-target (k*, n_real, R_p) under `key`, using debench's own walk."""
    out = []
    for p in range(len(T.perts)):
        pool = np.flatnonzero(T.S.real_sig[p])
        if pool.size == 0:
            continue
        if key is None:
            sig = T.S.pred_sig[p][pool]
            order = np.lexsort((T.genes[pool], -np.abs(T.S.pred_lfc[p, pool]),
                                T.S.pred_p[p, pool], T.S.pred_padj[p, pool], ~sig))
        else:
            order = np.lexsort((T.genes[pool], -np.asarray(key[p])[pool]))
        pool = pool[order]
        pool = pool[T.S.adj[p, pool]]
        n_real = int(T.S.real_sig[p].sum())
        if pool.size == 0:
            out.append((0, n_real, 0.0))
            continue
        ok = np.sign(T.S.pred_lfc[p, pool]) == np.sign(T.S.real_lfc[p, pool])
        pur = np.cumsum(ok) / np.arange(1, pool.size + 1)
        deep = np.flatnonzero(pur >= 0.9)
        k = int(deep[-1] + 1) if deep.size else 0
        out.append((k, n_real, k / n_real))
    return np.asarray(out, dtype=np.float64)


def main() -> None:
    T = H1Test()
    print("known-answer rows:")
    check(T)
    hf = np.load(os.path.join(HERE, "work", "h1test-feats.npz"))
    keys = {"production": None, "|H|": np.abs(T.H),
            "|cons| un-whitened": np.abs(hf["cons"]).astype(np.float64),
            "ORACLE": T.key_oracle()}
    print(f"\n{'key':<20}{'reach':>9}{'k*=0':>7}{'k*=1-2':>8}{'k*=3-9':>8}{'k*>=10':>8}"
          f"{'median k*':>11}{'max k*':>8}{'top-3 targets share':>21}")
    for nm, k in keys.items():
        A = kstar(T, k)
        ks, R = A[:, 0], A[:, 2]
        share = np.sort(R)[::-1][:3].sum() / max(R.sum(), 1e-9)
        print(f"{nm:<20}{R.mean():>9.5f}{int((ks == 0).sum()):>7}"
              f"{int(((ks >= 1) & (ks <= 2)).sum()):>8}{int(((ks >= 3) & (ks <= 9)).sum()):>8}"
              f"{int((ks >= 10).sum()):>8}{np.median(ks):>11.1f}{int(ks.max()):>8}"
              f"{share:>20.1%}")
    A = kstar(T, np.abs(T.H))
    n = A[:, 1]
    print("\nwhere the mass sits, under |H| — targets bucketed by n_real:")
    print(f"{'n_real':<14}{'targets':>9}{'mean k*':>10}{'mean R_p':>10}{'share of reach':>16}")
    tot = A[:, 2].sum()
    for lo, hi in ((1, 20), (21, 50), (51, 100), (101, 300), (301, 1000), (1001, 10 ** 9)):
        m = (n >= lo) & (n <= hi)
        if not m.any():
            continue
        print(f"{f'{lo}-{hi if hi < 10**9 else 'inf'}':<14}{int(m.sum()):>9}"
              f"{A[m, 0].mean():>10.2f}{A[m, 2].mean():>10.4f}"
              f"{A[m, 2].sum()/tot:>15.1%}")
    print("\nprefix accuracy needed for k* to reach a given depth (>=90% purity):")
    print("   k*=1 needs 1/1 correct · k*=5 needs 5/5 · k*=10 needs 9/10 · k*=20 needs 18/20")
    print("   at the measured front accuracy of ~0.62 per gene, P(5/5) ~ 0.09 if independent —")
    print("   so `reach` is a bet on a small number of targets whose front is nearly clean,")
    print("   which is exactly why a key that is uniformly decent loses to one that is")
    print("   excellent in its top few.")


if __name__ == "__main__":
    main()
