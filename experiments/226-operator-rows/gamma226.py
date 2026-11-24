"""226 Phase 0 gates KA-4/KA-5 — reproduce 193's vertex and 212's exponent result, AND
measure the gamma axis's own cosine/pds trade rate, which is this experiment's yardstick.

Why the yardstick matters more than the gates
---------------------------------------------
212 closed gamma by a RATIO argument, and said so explicitly: "a ratio is instrument-
independent".  Along the gamma axis the operator buys transfer cosine and `pds` pays for it at
`d(overall pds)/d(cos) = -0.3517` against `d(overall mse*)/d(cos) = +0.3655` -- they cancel to
3.8%, which is why gamma is a knife-edge and closed.

The row axis is a DIFFERENT axis.  If more fitting rows raise cosine at the same `pds` cost per
unit cosine as gamma does, then "more rows" is just gamma in disguise and 212 has already closed
it.  If the row axis trades at a materially better rate -- or does not charge `pds` at all --
then it is a new lever and has to be priced on its own.  So the number this script has to
produce is not a gate, it is `d(pds)/d(cos)` along gamma, measured on the SAME fold, the SAME
render and the SAME code path the row curve uses.

Gates (rule 10), all against 212's published `work/plolo212.json` rows on `kolf2`:
  KA-4  gamma 0.65, a15 render:  pds 0.611201  ==  193 phase2's shipped-operator row
  KA-5  gamma 0.65, a15 render:  cos 0.0078207 ==  212's `a15_cos`
  KA-6  cosine is monotone DECREASING in gamma over the whole grid  (212 §5)
  KA-7  a23 `pds` is monotone INCREASING in gamma from 0 to 0.80    (212 §5)
"""
from __future__ import annotations

import gc
import json
import os
import resource
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "src"))

import plolo                                                       # noqa: E402
from plolo import build_fold, endpoints                            # noqa: E402
from run_arms import CONFIGS                                       # noqa: E402
from rows226 import (PROD_OP, PROD_SHA, cd4_cached, frames_for,    # noqa: E402
                     score, sha256_file)

GRID = [0.00, 0.30, 0.35, 0.40, 0.45, 0.50, 0.525, 0.55, 0.575, 0.60, 0.625,
        0.65, 0.70, 0.75, 0.80]
EXPECT = {  # 212 work/plolo212.json, fold kolf2
    0.65: dict(a15_pds=0.611201, a15_cos=0.0078207, a23_pds=0.610128, a23_cos=0.0075590),
    0.30: dict(a15_pds=0.608216, a15_cos=0.0085890, a23_pds=0.602443, a23_cos=0.0081191),
    0.00: dict(a15_pds=0.588536, a15_cos=0.0092404, a23_pds=0.579809, a23_cos=0.0084958),
}


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30


def main() -> None:
    t0 = time.time()
    got = sha256_file(PROD_OP)
    assert got == PROD_SHA, f"rule 11: operator sha {got[:16]} != {PROD_SHA[:16]}"
    print(f"operator {PROD_OP.relative_to(ROOT)}  sha {got[:16]}  ASSERTED", flush=True)

    fold = build_fold("kolf2")
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    wcols = np.asarray(fold["wcols"])
    cd4 = cd4_cached(fold["targets"])
    frames = frames_for(fold)

    with np.load(PROD_OP) as z:
        evals = np.asarray(z["evals"], dtype=np.float64)
        evecs = np.asarray(z["evecs"], dtype=np.float64)
    print(f"loaded eigensystem {evecs.shape}, RSS peak {peak_gib():.2f} GiB", flush=True)

    rows = []
    for g in GRID:
        for arm in ("a15", "a23"):
            cfg = dict(CONFIGS[arm])
            cfg["gamma"] = g
            r = score(fold, cfg, evals, evecs, cd4, frames, truth, wcols)
            rows.append(dict(gamma=g, arm=arm, pds=r["chal_pds"], cos=r["chal_cos"],
                             cos_agg=r["chal_cos_agg"], per_frame={
                                 k: {kk: vv for kk, vv in v.items()
                                     if not kk.startswith(("cos_per", "pds_per", "cos_wrong"))}
                                 for k, v in r["per_frame"].items()}))
        a15 = rows[-2]
        a23 = rows[-1]
        print(f"[g={g:.3f}] a15 pds {a15['pds']:.6f} cos {a15['cos']:.7f} | "
              f"a23 pds {a23['pds']:.6f} cos {a23['cos']:.7f}", flush=True)
    del evals, evecs
    gc.collect()

    # ---- gates -----------------------------------------------------------------------
    gates = {}
    for g, exp in EXPECT.items():
        for arm in ("a15", "a23"):
            r = [x for x in rows if x["gamma"] == g and x["arm"] == arm][0]
            for key, field in (("pds", "pds"), ("cos", "cos")):
                want = exp[f"{arm}_{key}"]
                dev = abs(r[field] - want)
                gates[f"g{g:g}_{arm}_{key}"] = dict(expected=want, measured=r[field],
                                                    abs_dev=dev, passed=bool(dev < 5e-6))
    a15c = [x["cos"] for x in rows if x["arm"] == "a15"]
    a23p = [x["pds"] for x in rows if x["arm"] == "a23"]
    gates["KA6_cos_monotone_decreasing"] = dict(
        passed=bool(all(b < a for a, b in zip(a15c, a15c[1:]))),
        first=a15c[0], last=a15c[-1], ratio=a15c[0] / a15c[-1])
    # 212 §5 claims `pds` is maximised at gamma = 0.80 "monotone the other way" across its
    # FOUR folds.  On `kolf2` alone its own published rows are not strictly monotone -- they
    # dip at gamma = 0.70 (0.610128 -> 0.610088) -- so the single-fold gate is the endpoint
    # claim plus a count of inversions, and both are checked against 212's own numbers.
    inversions = sum(1 for a, b in zip(a23p, a23p[1:]) if b <= a)
    gates["KA7_a23_pds_max_at_080_on_kolf2"] = dict(
        passed=bool(a23p[-1] == max(a23p) and inversions <= 1),
        first=a23p[0], last=a23p[-1], argmax_gamma=GRID[int(np.argmax(a23p))],
        inversions=inversions,
        note="212's own kolf2 rows carry exactly one inversion, at gamma 0.70")

    # ---- the yardstick: the gamma axis's own trade rate, local to production ----------
    def at(g, arm, f):
        return [x for x in rows if x["gamma"] == g and x["arm"] == arm][0][f]

    trade = {}
    for arm in ("a15", "a23"):
        for lo, hi in ((0.625, 0.65), (0.60, 0.70), (0.50, 0.80)):
            dcos = at(hi, arm, "cos") - at(lo, arm, "cos")
            dpds = at(hi, arm, "pds") - at(lo, arm, "pds")
            trade[f"{arm}_{lo:g}_{hi:g}"] = dict(
                d_cos=dcos, d_pds=dpds,
                d_pds_per_d_cos=(dpds / dcos if dcos else None),
                rel_d_cos_pct=100 * dcos / at(lo, arm, "cos"))
    out = dict(operator=str(PROD_OP.relative_to(ROOT)), operator_sha256=got,
               fold="kolf2", grid=GRID, rows=rows, gates=gates, gamma_trade=trade,
               peak_rss_gib=round(peak_gib(), 2), seconds=round(time.time() - t0, 1))
    (HERE / "work/gamma226.json").write_text(json.dumps(out, indent=1))

    print("\n=== gates ===")
    for k, v in gates.items():
        print(f"  {'PASS' if v['passed'] else 'FAIL'}  {k}" +
              (f"  exp {v['expected']:.7f} got {v['measured']:.7f} |d| {v['abs_dev']:.2e}"
               if "expected" in v else ""))
    print("\n=== the gamma axis's trade rate (this experiment's yardstick) ===")
    for k, v in trade.items():
        print(f"  {k:18s} dcos {v['d_cos']:+.3e} ({v['rel_d_cos_pct']:+.3f}%)  "
              f"dpds {v['d_pds']:+.3e}  dpds/dcos {v['d_pds_per_d_cos']:+.4f}")
    print(f"\nwrote work/gamma226.json ({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)")
    if not all(v["passed"] for v in gates.values()):
        raise SystemExit("PHASE 0 GATE FAILED")


if __name__ == "__main__":
    main()
