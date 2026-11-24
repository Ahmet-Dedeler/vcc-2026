"""196 — `pds` under emitted-depth quantisation, on PLOLO (162's instrument).

162 classifies `pds` as the one member that is locally priceable, through PLOLO in the
challenge frame (LOO sigma 0.00149 raw, sign 6/6).  188 Phase 3 forbids sizing a `pds` arm
with `debench.Surface.pds()`; this file uses PLOLO and nothing else.

The arm is a pure quantisation arm: 188's `deliver()` renderer is kept VERBATIM (mass-fix
included, so the known-answer row below is comparable), and the ONLY thing added is the step
production actually performs and 188 omits --

    T = floor(n L p),   q = T / sum(T),   delivered = log1p(TS q) - b_ctrl

which is the emitter's integer rounding.  `m = inf` is 188's own row, i.e. the ceiling an
infinitely deep emitter would reach, and `m = 1` is production's realised depth.

KNOWN-ANSWER ROW (rule 10, abort on failure): at `m = inf` the four folds must reproduce
169's cached PLOLO `chalA/delivered` a23 rows -- 0.790712, 0.790896, 0.771621, 0.615317.

    uv run python pds196.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
E188 = ROOT / "experiments/188-magnitude-calibration"
sys.path.insert(0, str(E188))
sys.path.insert(0, str(ROOT / "experiments/169-sharpening-six-member"))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
from budget188 import FOLDS, PDS_SPAN, TS, deliver, load, mass_consistent_effect, pds_of  # noqa: E402

WORK = HERE / "work"
NL_1X = 400 * 20_109            # context A's n * L, production's realised depth
DEPTHS = tuple(float(x) for x in os.environ.get("DEPTHS", "1,2,4,8,16,64").split(","))
KA = {"k562_gwps": 0.790712, "orion_hct116": 0.790896,
      "orion_hek293t": 0.771621, "kolf2": 0.615317}


def deliver_q(effects, b, gate, nL):
    """188's `deliver`, plus the emitter's integer rounding at total depth `nL`."""
    out = np.empty_like(effects)
    for i in range(effects.shape[0]):
        g = mass_consistent_effect(np.where(gate, effects[i], 0.0), b, gate, "minnorm")
        p = np.clip(np.expm1(b + g), 0.0, None) / TS
        p = p / float(p.sum())
        T = np.floor(nL * p)
        out[i] = np.log1p(TS * (T / T.sum())) - b
    return out


def main() -> None:
    z = np.load(E188 / "work/chalframes.npz")
    bA, gA = z["bA"], z["gA"]
    F = {f: load(f) for f in FOLDS}
    res, quant = {}, {}
    print("PLOLO, challenge frame A, delivered — `pds` against emitted depth\n")
    for f in FOLDS:
        P0 = np.where(gA[None, :], F[f]["pred"], 0.0)
        Tsc = F[f]["truth"][:, F[f]["scored"]]
        inf = deliver(P0, bA, gA)
        res[(f, "inf")] = pds_of(inf[:, F[f]["scored"]], Tsc)
        for m in DEPTHS:
            D = deliver_q(P0, bA, gA, NL_1X * m)
            res[(f, m)] = pds_of(D[:, F[f]["scored"]], Tsc)
            # how much of the delivered effect's ENERGY is rounding noise
            num = float(np.linalg.norm(D[:, F[f]["scored"]] - inf[:, F[f]["scored"]]))
            den = float(np.linalg.norm(inf[:, F[f]["scored"]]))
            quant[(f, m)] = num / den
        print(f"  [{f}] {F[f]['pred'].shape[0]} perturbations, "
              f"{F[f]['scored'].size} scored genes  done", flush=True)

    print("\n--- KNOWN-ANSWER GATE: m = inf must reproduce 169's cached PLOLO chalA rows ---")
    worst = 0.0
    for f in FOLDS:
        w = KA[f]
        g = res[(f, "inf")]
        worst = max(worst, abs(g - w))
        print(f"  {f:<16} 169 {w:.6f}   here {g:.6f}   |d| {abs(g-w):.2e}")
    print(f"  worst |d pds| = {worst:.2e}")
    if worst > 2e-6:
        raise SystemExit("ABORT (rule 10): PLOLO identity row does not reproduce")

    print(f"\n{'depth':>7}{'n*L':>12}{'counts/cell':>13}" +
          "".join(f"{f[:11]:>13}" for f in FOLDS) + f"{'MEAN pds':>11}"
          f"{'d vs 1x':>10}{'scaled/6':>10}{'||round||/||eff||':>19}")
    rows = []
    base = float(np.mean([res[(f, 1.0)] for f in FOLDS]))
    for m in list(DEPTHS) + ["inf"]:
        vals = [res[(f, m)] for f in FOLDS]
        mu = float(np.mean(vals))
        qn = (float(np.mean([quant[(f, m)] for f in FOLDS])) if m != "inf" else 0.0)
        d = mu - base
        rows.append(dict(depth=str(m), pds=mu, d_pds=d, scaled_over6=d / PDS_SPAN / 6,
                         quant_energy=qn, per_fold={f: res[(f, m)] for f in FOLDS}))
        lab = "inf" if m == "inf" else f"{m:g}"
        nl = "-" if m == "inf" else f"{NL_1X*m:,.0f}"
        cc = "-" if m == "inf" else f"{20109*m:,.0f}"
        print(f"{lab:>7}{nl:>12}{cc:>13}" + "".join(f"{v:>13.6f}" for v in vals) +
              f"{mu:>11.6f}{d:>+10.6f}{d/PDS_SPAN/6:>+10.6f}{qn:>19.5f}")
    print(f"\n  span used for scaling: PDS_SPAN = {PDS_SPAN:.5f}; 162's PLOLO LOO sigma is "
          f"0.00149 RAW pds, i.e. {0.00149/PDS_SPAN/6:.6f} overall.")
    json.dump(rows, open(WORK / "pds196.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
