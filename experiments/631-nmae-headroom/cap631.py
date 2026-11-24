#!/usr/bin/env python3
"""631 cap — can a DOWN-to-ballast pass deliver the head-DOWN magnitude lever (m = 2) that decomp631 found?

decomp631: called HEAD genes delivered DOWN carry 12-16 % of H1's nmae numerator; their signs are right (perfect signs
buy -0.007) but their magnitudes are too small (|T| oracle -0.11 / -0.15 raw); one cross-fitted global factor m = 2
buys -0.017 (f 1) / -0.022 (f 0.35) raw. The emitter parks a DOWN gene's total on mass-holder carriers above the
control max, which inflates the arithmetic mean. Lowering the mean at a fixed total means moving counts carrier ->
ballast (a count is worth cw[a] - cw[b] ~ 1/(n s_a) of mean CPM), the mirror of 630.
Per zs-applied seed-0 H1 block and each head gene with delivered lfc < 0 (lib521's design head), this measures:
  need     counts to reach lfc = 2 x today's (ΔA = A - A_target at the median carrier-to-ballast gain);
  carrier  counts the carriers can give with ZERO rank change: above-cmax holders down to the minimum above cmax, plus
           in-interval room (cell_bounds x_lo) on the rest;
  ballast  counts the ballast rows can take: zero-cost room (x_hi, open above the max) and, separately, with the gene's
           net U budget (|ΔAUC| <= 1e-4, upward crossings only) spent;
  currency the ballast rows' sub-gate counts (they pay the carriers' library back; the mirror of 630's currency).
Output: per target and mean, the fraction of the m = 2 need reachable (zero-cost / budgeted), and need vs currency.
Analysis only, ~4 GB.   uv run python experiments/631-nmae-headroom/cap631.py [<n_targets>]
"""
from __future__ import annotations

import glob
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/630-rank-neutral-up"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import rnup630 as U  # noqa: E402
import src_rank_up as M  # noqa: E402
import lib521 as L  # noqa: E402


def max_rise_k(x, s_b, c, budget):
    """Largest k >= 0 with ranks_of(x + k) - ranks_of(x) <= budget (x2 units); open above the control max."""
    r0 = int(U.ZS.ranks_of(np.array([x]), np.array([s_b]), c)[0])
    top = 2 * c.size
    if r0 >= top:                       # already above every control: unlimited
        return 1 << 30
    lo, hi = 0, 1
    while int(U.ZS.ranks_of(np.array([x + hi]), np.array([s_b]), c)[0]) - r0 <= budget and hi < (1 << 30):
        lo, hi = hi, hi * 2
        if int(U.ZS.ranks_of(np.array([x + lo]), np.array([s_b]), c)[0]) >= top:
            return 1 << 30
    while lo < hi - 1:
        mid = (lo + hi) // 2
        if int(U.ZS.ranks_of(np.array([x + mid]), np.array([s_b]), c)[0]) - r0 <= budget:
            lo = mid
        else:
            hi = mid
    return lo


def main() -> int:
    n_t = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    dauc = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-4   # per-gene |ΔAUC| budget for ballast crossings
    scope = sys.argv[3] if len(sys.argv) > 3 else "head"          # head | called (all called DOWN genes)
    r, c, _ = U.E.load_reference()
    keep = np.asarray(c.keep, dtype=bool)
    tie = M.scorer_tie(c)
    S = L.load_state(0)
    cs, ref = c.cs, np.asarray(c.ref_lin, float)
    files = sorted(glob.glob(str(ROOT / f"experiments/630-rank-neutral-up/work/blocks/{U.TAG.format(seed=0)}_*.npz")))[:n_t]
    zspec = U.ZS.parse_spec("dn_k250_zs")
    tot = []
    for i, f in enumerate(files):
        X0, n_car, target = U.load_block(f)
        ti = list(S["targets"]).index(target)
        _, _, _, _, l0 = U.E.block_stats(X0, c)
        Xz, _ = U.ZS.postpass(X0, n_car, c, l0, zspec)
        p, q, l = M.de_stats(Xz, c, tie)
        n = Xz.shape[0]
        s = Xz.sum(1) / 1e6
        cw = 1.0 / (n * s)
        car, bal = np.arange(n_car), np.arange(n_car, n)
        Xk = Xz[:, keep]
        A = (Xk / s[:, None]).mean(0)
        sel0 = S["head"][ti] if scope == "head" else np.ones(l.shape, bool)
        heads = np.flatnonzero(sel0 & np.isfinite(l) & (l < 0) & (q < 0.05))
        gain = float(np.median(cw[car]) - np.median(cw[bal]))
        capU = int(np.floor(dauc * 2.0 * n * cs.shape[1]))
        need = car_side = bal0 = balU = 0
        reach0 = reachU = 0.0
        for j in heads:
            A_t = (ref[j] + 1e-9) * 2.0 ** (2.0 * l[j]) - 1e-9           # m = 2: lfc -> 2 x today's
            nd = int(np.ceil(max(A[j] - A_t, 0.0) / gain))
            cmaxj = cs[j][-1]
            xc = Xk[car, j]
            above = np.log1p(xc / s[car]) > cmaxj
            xl, *_ = U.ZS.cell_bounds(xc, s[car], cs[j])
            free = 0
            for a in np.flatnonzero(above):
                free += int(xc[a]) - M._min_above(cmaxj, s[a])
            free += int(np.maximum(xc[~above] - xl[~above], 0).sum())
            _, xh, *_ = U.ZS.cell_bounds(Xk[bal, j], s[bal], cs[j])
            z = int(np.where(xh >= U.ZS.BIG, 1 << 30, np.maximum(xh - Xk[bal, j], 0)).sum())
            budget, bu = capU, 0
            for b in np.argsort(-Xk[bal, j]):
                k = max_rise_k(int(Xk[bal[b], j]), float(s[bal[b]]), cs[j], budget)
                d = int(U.ZS.ranks_of(np.array([Xk[bal[b], j] + min(k, 1 << 20)]), np.array([s[bal[b]]]), cs[j])[0]) \
                    - int(U.ZS.ranks_of(np.array([Xk[bal[b], j]]), np.array([s[bal[b]]]), cs[j])[0])
                bu += k
                budget -= d
                if budget <= 0:
                    break
            need += nd
            car_side += free
            bal0 += min(z, 1 << 30)
            balU += min(bu, 1 << 30)
            if nd:
                reach0 += min(1.0, min(free, z) / nd)
                reachU += min(1.0, min(free, bu) / nd)
        sub_bal = int(Xz[np.ix_(bal, np.flatnonzero(~keep))].sum())
        k = max(len(heads), 1)
        rec = dict(target=target, head_down=len(heads), need=need, car=car_side, bal0=min(bal0, 10**9),
                   balU=min(balU, 10**9), frac0=reach0 / k, fracU=reachU / k, sub_bal=sub_bal, n_bal=int(bal.size))
        tot.append(rec)
        print(f"{target:10s} head-DOWN called {rec['head_down']:3d} | need {need:7d} | carrier {car_side:7d} | ballast "
              f"zero-cost {rec['bal0']:>10} budgeted {rec['balU']:>10} | reachable frac zero-cost {rec['frac0']:.2f} "
              f"budgeted {rec['fracU']:.2f} | ballast sub-gate {sub_bal} (K {bal.size})", flush=True)
    m = {k: float(np.mean([t[k] for t in tot])) for k in tot[0] if k != "target"}
    print("\nMEAN " + " ".join(f"{k} {v:.4g}" for k, v in m.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
