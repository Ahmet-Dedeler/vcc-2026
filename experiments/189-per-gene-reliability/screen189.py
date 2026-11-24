"""189 Phase 2.1 — the prior as an ORDERING KEY on h1-test, against `|H|`'s 0.23533.

`reach` reads only the FRONT of the walk (172 §4: rank correlation over the pool is
anti-correlated with `reach` on this surface), so top-1% / top-5% pooled and macro sign
accuracy are reported next to every `reach`.

Every key here is `|H|` re-weighted by a GENE-LEVEL scalar, i.e. `log key = log|H| + lam*theta
+ mu*z(expression)`.  A gene-level term is constant within a target, so it cannot change the
BETWEEN-target base rates the reach walk averages -- only the within-target order, which is the
only thing `reach` lets us move (171 §0.1).

Selection honesty: `lam` and `mu` are swept ON h1-test.  Any positive row is therefore an upper
bound selected on the evaluation surface, and is read against 172 §4.1's paired-bootstrap SE
for a reach difference (0.005-0.016).  The permutation null and the expression-only control are
run on the same lines of code.
"""
from __future__ import annotations

import json
import os

import numpy as np

from common189 import WORK
from harness189 import H_ORACLE, Surface189, audit, check, paired_bootstrap

LAMS = [-4, -2, -1, -0.5, -0.25, 0.25, 0.5, 1, 2, 4]
MUS = [-1.0, -0.5, -0.25, 0.25, 0.5, 1.0]


def zs(x: np.ndarray) -> np.ndarray:
    m = np.isfinite(x)
    o = np.zeros_like(x, dtype=np.float64)
    o[m] = (x[m] - x[m].mean()) / max(x[m].std(), 1e-12)
    return o


def main() -> None:
    T = Surface189()
    check(T)
    print()
    audit(T)

    absH = np.abs(T.H)
    th = T.theta
    ze = zs(np.log1p(np.expm1(T.b_ctrl)))     # monotone in expression; z-scored for the sweep
    base = np.log(np.maximum(absH, 1e-30))

    keys: dict[str, np.ndarray] = {
        "production order": None,
        "|H| (INCUMBENT)": absH,
        "theta alone (gene-only)": np.broadcast_to(th, absH.shape).copy(),
        "expression alone (gene-only)": np.broadcast_to(ze, absH.shape).copy(),
        "|H| x raw r_g": absH * np.nan_to_num(T.raw, nan=float(np.nanmean(T.raw))),
    }
    for lam in LAMS:
        keys[f"|H| x exp({lam:+.2f} theta)"] = np.exp(base + lam * th)
    for mu in MUS:
        keys[f"|H| x exp({mu:+.2f} z_expr)"] = np.exp(base + mu * ze)

    ref = T.reach_vec(keys["|H| (INCUMBENT)"])
    rows = []
    print(f"\n{'=' * 118}\nORDERING KEYS on h1-test  "
          f"(incumbent |H| 0.23533, H-oracle {H_ORACLE:.5f}, full oracle 0.64553)")
    print(f"{'key':<34}{'reach':>9}{'d vs |H|':>10}{'SE':>8}{'P(bet)':>8}"
          f"{'top1%':>8}{'top5%':>8}{'top10%':>8}{'macro1%':>9}{'d1':>7}{'d10':>7}")
    for name, k in keys.items():
        r = T.reach(k) if k is not None else T.reach()
        v = T.reach_vec(k) if k is not None else T.reach_vec()
        d, se, lo, hi, pb = paired_bootstrap(v, ref)
        tx = T.topx(k) if k is not None else None
        mx = T.topx_macro(k) if k is not None else None
        dc = T.deciles(k) if k is not None else None
        rows.append(dict(key=name, reach=r, d=d, se=se, ci=[lo, hi], p_better=pb,
                         topx=tx, macro=mx, deciles=dc))
        if tx is None:
            print(f"{name:<34}{r:>9.5f}{d:>+10.5f}{se:>8.5f}{pb:>8.3f}"
                  f"{'':>8}{'':>8}{'':>8}{'':>9}{'':>7}{'':>7}")
        else:
            print(f"{name:<34}{r:>9.5f}{d:>+10.5f}{se:>8.5f}{pb:>8.3f}"
                  f"{tx[0]:>8.3f}{tx[2]:>8.3f}{tx[3]:>8.3f}{mx[0]:>9.3f}"
                  f"{dc[0]:>7.3f}{dc[9]:>7.3f}")

    # ---- the (lam, mu) grid ----------------------------------------------------------
    print(f"\n{'=' * 118}\nJOINT GRID  reach for  |H| x exp(lam*theta + mu*z_expr)")
    print(f"{'lam \\ mu':>10}" + "".join(f"{m:>10.2f}" for m in [0.0] + MUS))
    grid = {}
    for lam in [0.0] + LAMS:
        line = f"{lam:>10.2f}"
        for mu in [0.0] + MUS:
            r = T.reach(np.exp(base + lam * th + mu * ze))
            grid[f"{lam}|{mu}"] = r
            line += f"{r:>10.5f}"
        print(line)

    # ---- permutation null -------------------------------------------------------------
    print(f"\n{'=' * 118}\nPERMUTATION NULL — theta shuffled ACROSS GENES (10 seeds)")
    best_lam = max(LAMS, key=lambda l: T.reach(np.exp(base + l * th)))
    print(f"h1-selected lam = {best_lam:+.2f}, real theta reach "
          f"{T.reach(np.exp(base + best_lam * th)):.5f}")
    null = []
    for s in range(10):
        rng = np.random.default_rng(1000 + s)
        null.append(T.reach(np.exp(base + best_lam * rng.permutation(th))))
    print(f"null over 10 seeds: mean {np.mean(null):.5f} sd {np.std(null, ddof=1):.5f} "
          f"range [{min(null):.5f}, {max(null):.5f}]")
    print(f"real minus null mean: {T.reach(np.exp(base + best_lam*th)) - np.mean(null):+.5f} "
          f"({(T.reach(np.exp(base + best_lam*th)) - np.mean(null))/np.std(null, ddof=1):+.2f} "
          f"null sd)")

    json.dump(dict(rows=rows, grid=grid, null=null, best_lam=best_lam,
                   incumbent=0.23533, h_oracle=H_ORACLE, oracle=0.64553),
              open(os.path.join(WORK, "screen189.json"), "w"), indent=2)
    print("\nwrote work/screen189.json")


if __name__ == "__main__":
    main()
