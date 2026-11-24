"""202 — read the EMITTED h1-test tables with 192's own instruments.

`ka202_prescribed.py` reproduced 192 §5.3's grid exactly, but that grid PRESCRIBES a call
set on 185's frozen tables.  This scores the tables the patched emitter actually produced,
so the four DE members are read off the emitted order and the emitted call set -- the
`_order(key=None)` path, i.e. "significant calls first, then ascending p_adj, then
ascending p, then descending |lfc|, then gene name", which is what the scorer does.

Both arms are emitted with production's single-stream seeding, so the pair is paired: the
same seeds, the same control block, the same per-gene totals (digest-identical on 85/85),
one variable.

`n_real` is scaled PROPORTIONALLY (`PowerSurface`), never `Surface(n_real_target=N)` --
harness trap 2.  Reported at 192's own four levels.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
sys.path.insert(0, f"{ROOT}/experiments/192-predicted-callset")
from debench import Surface  # noqa: E402
from frontier192 import INV, PowerSurface  # noqa: E402

ARMS = ("rc_hybrid", "rc_callset")
# 192's own prescribed-set prediction for this arm, f = 0.50 (`ka202_prescribed.log`)
PRESCRIBED = dict(fid=+0.00097, jac=+0.01288, reach=+0.01123)


def main() -> None:
    paths = {a: os.path.join(WORK, f"de202.{a}.npz") for a in ARMS}
    base = Surface(paths["rc_hybrid"])
    nr = np.array([int(base.real_sig[p].sum()) for p in range(len(base.perts))])
    print(f"h1-test emitted tables: {len(base.perts)} perturbations x {base.genes.size} "
          f"gated genes; reference n_real median {np.median(nr):.0f} mean {nr.mean():.0f}")

    out = {}
    for f in (0.35, 0.40, 0.50, 0.65):
        row = {}
        for arm in ARMS:
            s = PowerSurface(paths[arm], f)
            npd = np.array([int((s.pred_sig[p] & s.adj[p]).sum())
                            for p in range(len(s.perts))])
            row[arm] = dict(fid=s.fid(), jac=s.jac(), reach=s.reach(), nmae=s.nmae(),
                            n_pred=float(np.median(npd)))
        d = {m: row["rc_callset"][m] - row["rc_hybrid"][m]
             for m in ("fid", "jac", "reach", "nmae")}
        overall = ((d["fid"] * INV["fid"] + d["jac"] * INV["jac"]) / 6
                   + d["reach"] * INV["reach"] / 6 * 1.290)
        overall_c = ((d["fid"] * INV["fid"] + d["jac"] * INV["jac"]) / 6
                     + d["reach"] * INV["reach"] / 6 * 0.4277)
        out[f"{f}"] = dict(rows=row, delta=d, sum_129=overall, sum_043=overall_c)
        nrf = np.array([int(PowerSurface(paths["rc_hybrid"], f).real_sig[p].sum())
                        for p in range(len(base.perts))])
        print(f"\nf = {f}  (median n_real {np.median(nrf):.0f}, mean {nrf.mean():.0f})")
        print(f"  {'arm':<12}{'n_pred':>8}{'fid':>10}{'jac':>10}{'reach':>10}{'nmae':>10}")
        for arm in ARMS:
            r = row[arm]
            print(f"  {arm:<12}{r['n_pred']:>8.0f}{r['fid']:>10.5f}{r['jac']:>10.5f}"
                  f"{r['reach']:>10.5f}{r['nmae']:>10.5f}")
        print(f"  {'DELTA':<12}{'':>8}{d['fid']:>+10.5f}{d['jac']:>+10.5f}"
              f"{d['reach']:>+10.5f}{d['nmae']:>+10.5f}"
              f"   overall @1.290 {overall:+.5f}  @0.4277 {overall_c:+.5f}")

    print("\nEMITTED vs 192's PRESCRIBED prediction, f = 0.50 (raw h1-test)")
    d = out["0.5"]["delta"]
    for m in ("fid", "jac", "reach"):
        print(f"  {m:<6} emitted {d[m]:+.5f}   192 prescribed {PRESCRIBED[m]:+.5f}   "
              f"ratio {d[m] / PRESCRIBED[m]:.2f}x")

    with open(os.path.join(WORK, "score202.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"\nwrote {WORK}/score202.json")


if __name__ == "__main__":
    main()
