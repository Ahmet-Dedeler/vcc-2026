"""189 — is there a per-GENE reliability channel on h1-test at all, and does theta point at it?

`members189.py` measured a gene-level ORACLE worth +0.163 `q_w` while the donor-fitted prior
moves it by +0.0003.  Two very different explanations:

  (i)  the oracle is real and our prior simply points the wrong way, or
  (ii) the oracle is mostly NOISE -- each gene appears in only a handful of the 85 targets'
       DE pools, so `acc_g` fitted on the same rows it is scored on is guaranteed to look good.

The falsifier is a SPLIT-HALF: estimate the per-gene accuracy on one half of the targets and
read it on the other half.  A gene-level channel that survives that is real and is the ceiling
185/188 should be told about; one that does not means the whole family is closed on this
surface and our prior's failure is not a failure of THIS prior.

Also reported: the direct correlation between the donor-fitted `theta_g` and h1-test's own
per-gene accuracy, which is the single number that says whether cross-line donor agreement and
prediction reliability are the same quantity.
"""
from __future__ import annotations

import json
import os

import numpy as np

from common189 import WORK, spearman
from harness189 import Surface189, check


def gene_acc(T, targets):
    """(k_g, n_g) over the reference DE pool, restricted to `targets`."""
    k = np.zeros(len(T.genes))
    n = np.zeros(len(T.genes))
    for p in targets:
        idx = np.flatnonzero(T.pool[p])
        k[idx] += T.correct[p, idx]
        n[idx] += 1
    return k, n


def main() -> None:
    T = Surface189()
    check(T)
    P = len(T.perts)
    k_all, n_all = gene_acc(T, range(P))
    acc = np.where(n_all > 0, k_all / np.maximum(n_all, 1), np.nan)
    live = n_all > 0
    print(f"\n{'=' * 96}\nGENE-LEVEL CHANNEL ON h1-test")
    print(f"  genes appearing in the DE pool at all: {int(live.sum())} of {len(T.genes)}")
    q = np.percentile(n_all[live], [50, 75, 90, 99])
    print(f"  targets per gene in the pool: median {q[0]:.0f}  p75 {q[1]:.0f}  "
          f"p90 {q[2]:.0f}  p99 {q[3]:.0f}  max {n_all.max():.0f}")
    print(f"  per-gene accuracy: mean {np.nanmean(acc):.4f} sd {np.nanstd(acc):.4f}")

    m = live & np.isfinite(T.raw)
    print(f"\n  spearman(theta_g, h1 per-gene accuracy)      {spearman(T.theta[m], acc[m]):+.4f}"
          f"   over {int(m.sum())} genes")
    print(f"  spearman(raw r_g, h1 per-gene accuracy)      {spearman(T.raw[m], acc[m]):+.4f}")
    print(f"  spearman(expression, h1 per-gene accuracy)   "
          f"{spearman(T.b_ctrl[m], acc[m]):+.4f}")
    hw = live & (n_all >= 10)
    print(f"  the same three, genes with >=10 pool targets ({int(hw.sum())} genes): "
          f"{spearman(T.theta[hw], acc[hw]):+.4f} / {spearman(T.raw[hw], acc[hw]):+.4f} / "
          f"{spearman(T.b_ctrl[hw], acc[hw]):+.4f}")

    # ---- split-half -------------------------------------------------------------------
    print(f"\n{'=' * 96}\nSPLIT-HALF: is the gene-level channel real, or is the ORACLE noise?")
    print(f"{'seed':>5}{'rho(acc_A, acc_B)':>20}{'d10 spread A->B':>18}"
          f"{'q_w oracle-in-fold':>21}{'q_w oracle-cross':>19}")
    rows = []
    q0_all = T.qw()
    K0 = T.mass_ratio()
    for s in range(12):
        rng = np.random.default_rng(500 + s)
        perm = rng.permutation(P)
        A, B = perm[: P // 2], perm[P // 2:]
        kA, nA = gene_acc(T, A)
        kB, nB = gene_acc(T, B)
        both = (nA >= 3) & (nB >= 3)
        aA = kA[both] / nA[both]
        aB = kB[both] / nB[both]
        rho = spearman(aA, aB)
        o = np.argsort(aA)
        d = max(1, len(o) // 10)
        spread = float(aB[o[-d:]].mean() - aB[o[:d]].mean())
        # q_w on fold B, weighting by an in-fold oracle vs a cross-fold (fold-A) oracle
        def qw_fold(fold, w):
            num = den = 0.0
            for p in fold:
                idx = np.flatnonzero(T.pool[p])
                if idx.size == 0:
                    continue
                wt = np.abs(T.S.pred_lfc[p, idx]) * w[idx]
                ok = (np.sign(T.S.pred_lfc[p, idx]) ==
                      np.sign(T.S.real_lfc[p, idx])).astype(float)
                num += float((wt * ok).sum())
                den += float(wt.sum())
            return num / max(den, 1e-300)
        wB = np.where(nB > 0, (kB / np.maximum(nB, 1)) > 0.5, True).astype(float) + 1e-6
        wA = np.where(nA > 0, (kA / np.maximum(nA, 1)) > 0.5, True).astype(float) + 1e-6
        rows.append(dict(seed=s, rho=rho, spread=spread,
                         qw_in=qw_fold(B, wB), qw_cross=qw_fold(B, wA),
                         qw_base=qw_fold(B, np.ones(len(T.genes)))))
        print(f"{s:>5}{rho:>+20.4f}{spread:>+18.4f}{rows[-1]['qw_in']:>21.4f}"
              f"{rows[-1]['qw_cross']:>19.4f}")
    mr = lambda k_: float(np.mean([r[k_] for r in rows]))
    print(f"{'MEAN':>5}{mr('rho'):>+20.4f}{mr('spread'):>+18.4f}{mr('qw_in'):>21.4f}"
          f"{mr('qw_cross'):>19.4f}   (unweighted baseline {mr('qw_base'):.4f})")
    dq = mr("qw_cross") - mr("qw_base")
    print(f"\n  HONEST gene-level q_w headroom (cross-fold oracle - baseline): {dq:+.5f}")
    print(f"    -> d nmae raw {-2 * K0 * dq:+.5f}, scaled {1.6532 * 2 * K0 * dq:+.5f}, "
          f"overall {1.6532 * 2 * K0 * dq / 6:+.6f}")
    print(f"  IN-FOLD oracle (what members189 reported): "
          f"{mr('qw_in') - mr('qw_base'):+.5f} -- inflated by "
          f"{(mr('qw_in') - mr('qw_base')) / max(dq, 1e-9):.1f}x")

    json.dump(dict(rows=rows, rho_theta=spearman(T.theta[m], acc[m]),
                   rho_raw=spearman(T.raw[m], acc[m]),
                   rho_expr=spearman(T.b_ctrl[m], acc[m]),
                   n_pool_targets=list(map(float, q)), q0=q0_all, K0=K0),
              open(os.path.join(WORK, "genelevel189.json"), "w"), indent=2)
    print("\nwrote work/genelevel189.json")


if __name__ == "__main__":
    main()
