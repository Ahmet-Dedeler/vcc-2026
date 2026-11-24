#!/usr/bin/env python3
"""630 — how big is a DOWN-gene compensation currency? (the capacity bound found by diagcap630.py)

The UP pass pays each count a recipient carrier receives with one count of a sub-gate gene, moved to the donor ballast row,
so libraries stay exact. Carriers hold only ~0.55 % of their library in sub-gate genes, and the pass spends 93 % of it.
Candidate new currency: a kept gene g whose delivered lfc is DOWN, on a carrier a that sits ABOVE g's control max (zs's
"mass holders"). Moving one count of g from a to a ballast row b:
  * keeps a above the control max if a stays >= the minimum count above it: zero U change on the carrier side;
  * raises b's value for g: zero U change if b stays inside its current control interval (cell_bounds' x_hi; BIG above
    the control max), otherwise it crosses controls and spends g's own U budget;
  * lowers g's mean CPM by (cw[a] - cw[b]) per count: MORE down magnitude, the direction zs is capacity-bound on.
Per block (zs-applied, seed-0 H1, first N targets) this reports:
  sub-gate pool in carriers (today's currency);
  DOWN currency, carrier side: sum over DOWN genes and above-cmax carriers of (count - minimum count above cmax);
  ballast side: each DOWN gene's zero-cost room summed over ballast rows (capped at the carrier side per gene);
  PAIRED zero-cost currency = sum over genes of min(carrier side, ballast room); split called / uncalled;
  the mean shift it would apply if all used: median |Δlfc| on the genes that carry it.
  uv run python experiments/630-rank-neutral-up/currency630.py [<n_targets>]
Analysis only, ~4 GB (the H1 reference).
"""
from __future__ import annotations

import glob
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rnup630 as U  # noqa: E402
import src_rank_up as M  # noqa: E402


def main() -> int:
    n_t = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    r, c, _ = U.E.load_reference()
    keep = np.asarray(c.keep, dtype=bool)
    tie = M.scorer_tie(c)
    cs = c.cs
    cmax = cs[:, -1]
    ref = np.asarray(c.ref_lin, float)
    files = sorted(glob.glob(str(HERE / f"work/blocks/{U.TAG.format(seed=0)}_*.npz")))[:n_t]
    zspec = U.ZS.parse_spec("dn_k250_zs")
    tot = []
    for f in files:
        X0, n_car, target = U.load_block(f)
        _, _, _, _, l0 = U.E.block_stats(X0, c)
        Xz, _ = U.ZS.postpass(X0, n_car, c, l0, zspec)
        p, q, l = M.de_stats(Xz, c, tie)
        n = Xz.shape[0]
        s = Xz.sum(1) / 1e6
        cw = 1.0 / (n * s)
        car, bal = np.arange(n_car), np.arange(n_car, n)
        Xk = Xz[:, keep]
        sub_pool = int(Xz[np.ix_(car, np.flatnonzero(~keep))].sum())
        down = np.isfinite(l) & (l < 0)
        V = np.log1p(Xk[car] / s[car, None])
        above = (V > cmax[None, :]) & down[None, :]
        # minimum count that keeps carrier a above cmax_g (vectorised _min_above with corrections)
        xm = np.floor(np.expm1(cmax)[None, :] * s[car, None]).astype(np.int64) + 1
        for _ in range(3):
            xm = xm - (np.log1p((xm - 1) / s[car, None]) > cmax[None, :])
        for _ in range(3):
            xm = xm + ~(np.log1p(xm / s[car, None]) > cmax[None, :])
        free_car = np.where(above, np.maximum(Xk[car] - xm, 0), 0).sum(0)          # per gene
        gd = np.flatnonzero(free_car > 0)
        room = np.zeros(Xk.shape[1], dtype=np.int64)
        for j in gd:
            _, xh, *_ = U.ZS.cell_bounds(Xk[bal, j], s[bal], cs[j])
            rj = np.where(xh >= U.ZS.BIG, free_car[j], xh - Xk[bal, j])
            room[j] = int(np.minimum(np.maximum(rj, 0), free_car[j]).sum())
        paired = np.minimum(free_car, room)
        called = q < 0.05
        # mean shift if every paired count moved from the cheapest carriers' cw to the ballast's cw (upper bound)
        a_g = (Xk / s[:, None]).mean(0)
        dA = paired * (np.median(cw[car]) - np.median(cw[bal]))
        dl = np.abs(np.log2(np.maximum(a_g - dA, 1e-9) + 1e-9) - np.log2(a_g + 1e-9))
        rec = dict(target=target, sub_pool=sub_pool, down_genes_with_holders=int(gd.size),
                   car_side=int(free_car.sum()), paired=int(paired.sum()), paired_called=int(paired[called].sum()),
                   paired_uncalled=int(paired[~called].sum()), genes_paired=int((paired > 0).sum()),
                   median_dlfc=float(np.median(dl[paired > 0])) if (paired > 0).any() else 0.0)
        tot.append(rec)
        print(f"{target:10s} sub-gate pool {sub_pool:6d} | DOWN genes with holders {rec['down_genes_with_holders']:5d}, "
              f"carrier side {rec['car_side']:7d}, PAIRED zero-cost {rec['paired']:7d} (called {rec['paired_called']}, "
              f"uncalled {rec['paired_uncalled']}) on {rec['genes_paired']} genes, median |Δlfc| if all used "
              f"{rec['median_dlfc']:.3f}", flush=True)
    m = {k: float(np.mean([t[k] for t in tot])) for k in tot[0] if k != "target"}
    print("\nMEAN " + " ".join(f"{k} {v:.4g}" for k, v in m.items()))
    print(f"PAIRED / sub-gate pool = {m['paired'] / m['sub_pool']:.2f}x")
    return 0


if __name__ == "__main__":
    sys.exit(main())
