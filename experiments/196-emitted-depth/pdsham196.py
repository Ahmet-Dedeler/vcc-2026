"""196 — `pds` for the ROUNDING RULE at production's depth, on PLOLO (162's instrument).

`pds196.py` swept the depth with the rule held at `floor`.  This sweeps the RULE with the
depth held at production's `n L`.  Same renderer, same folds, same known-answer gate: at
`m = inf` the four folds must reproduce 169's cached PLOLO chalA rows.

Three rules at the identical `n L`:

  floor   `T = floor(n L p)`                         -- what `emit_counts` ships
  ham     largest-remainder over `n L p`             -- floor, then the leftover counts to
                                                        the largest fractional remainders
  inf     no rounding at all                         -- the ceiling

188 Phase 3 forbids `debench.Surface.pds()` from sizing a `pds` arm; this file uses PLOLO
and nothing else.

    uv run python pdsham196.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
E188 = ROOT / "experiments/188-magnitude-calibration"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(E188))
sys.path.insert(0, str(ROOT / "experiments/169-sharpening-six-member"))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
from budget188 import FOLDS, PDS_SPAN, TS, deliver, load, mass_consistent_effect, pds_of  # noqa: E402
from ham196 import hamilton                                              # noqa: E402

WORK = HERE / "work"
NL_1X = 400 * 20_109
KA = {"k562_gwps": 0.790712, "orion_hct116": 0.790896,
      "orion_hek293t": 0.771621, "kolf2": 0.615317}


def deliver_r(effects, b, gate, nL, rule):
    out = np.empty_like(effects)
    err = np.zeros(effects.shape[0])
    for i in range(effects.shape[0]):
        g = mass_consistent_effect(np.where(gate, effects[i], 0.0), b, gate, "minnorm")
        p = np.clip(np.expm1(b + g), 0.0, None) / TS
        p = p / float(p.sum())
        x = nL * p
        T = np.floor(x) if rule == "floor" else hamilton(x).astype(np.float64)
        nz = x > 0
        err[i] = float(np.abs((T[nz] - x[nz]) / x[nz]).mean())
        out[i] = np.log1p(TS * (T / T.sum())) - b
    return out, float(err.mean())


def main() -> None:
    z = np.load(E188 / "work/chalframes.npz")
    bA, gA = z["bA"], z["gA"]
    F = {f: load(f) for f in FOLDS}
    res, quant, relerr = {}, {}, {}
    print("PLOLO, challenge frame A, delivered — `pds` against the ROUNDING RULE\n")
    for f in FOLDS:
        P0 = np.where(gA[None, :], F[f]["pred"], 0.0)
        Tsc = F[f]["truth"][:, F[f]["scored"]]
        inf = deliver(P0, bA, gA)
        res[(f, "inf")] = pds_of(inf[:, F[f]["scored"]], Tsc)
        for rule in ("floor", "ham"):
            D, e = deliver_r(P0, bA, gA, NL_1X, rule)
            res[(f, rule)] = pds_of(D[:, F[f]["scored"]], Tsc)
            relerr[(f, rule)] = e
            num = float(np.linalg.norm(D[:, F[f]["scored"]] - inf[:, F[f]["scored"]]))
            quant[(f, rule)] = num / float(np.linalg.norm(inf[:, F[f]["scored"]]))
        print(f"  [{f}] {F[f]['pred'].shape[0]} perturbations, "
              f"{F[f]['scored'].size} scored genes  done", flush=True)

    print("\n--- KNOWN-ANSWER GATE: m = inf must reproduce 169's cached PLOLO chalA rows ---")
    worst = 0.0
    for f in FOLDS:
        w, g = KA[f], res[(f, "inf")]
        worst = max(worst, abs(g - w))
        print(f"  {f:<16} 169 {w:.6f}   here {g:.6f}   |d| {abs(g-w):.2e}")
    print(f"  worst |d pds| = {worst:.2e}")
    if worst > 2e-6:
        raise SystemExit("ABORT (rule 10): PLOLO identity row does not reproduce")
    # and the floor row must reproduce pds196's own m=1 row
    prev = {r["depth"]: r for r in json.load(open(WORK / "pds196.json"))}
    for f in FOLDS:
        d = abs(res[(f, "floor")] - prev["1.0"]["per_fold"][f])
        print(f"  {f:<16} pds196 m=1 {prev['1.0']['per_fold'][f]:.6f}   here "
              f"{res[(f,'floor')]:.6f}   |d| {d:.2e}")
        if d > 2e-6:
            raise SystemExit("ABORT (rule 10): the floor row does not reproduce pds196")

    print(f"\n{'rule':>7}{'n*L':>12}" + "".join(f"{f[:11]:>13}" for f in FOLDS) +
          f"{'MEAN pds':>11}{'d vs floor':>12}{'scaled/6':>10}"
          f"{'||round||/||eff||':>19}{'mean |rel err|':>16}")
    rows = []
    base = float(np.mean([res[(f, "floor")] for f in FOLDS]))
    for rule in ("floor", "ham", "inf"):
        vals = [res[(f, rule)] for f in FOLDS]
        mu = float(np.mean(vals))
        qn = float(np.mean([quant[(f, rule)] for f in FOLDS])) if rule != "inf" else 0.0
        re = float(np.mean([relerr[(f, rule)] for f in FOLDS])) if rule != "inf" else 0.0
        d = mu - base
        rows.append(dict(rule=rule, pds=mu, d_pds=d, scaled_over6=d / PDS_SPAN / 6,
                         quant_energy=qn, rel_err=re,
                         per_fold={f: res[(f, rule)] for f in FOLDS}))
        nl = "-" if rule == "inf" else f"{NL_1X:,.0f}"
        print(f"{rule:>7}{nl:>12}" + "".join(f"{v:>13.6f}" for v in vals) +
              f"{mu:>11.6f}{d:>+12.6f}{d/PDS_SPAN/6:>+10.6f}{qn:>19.5f}{re:>16.6f}")
    print(f"\n  PDS_SPAN = {PDS_SPAN:.5f}; 162's PLOLO LOO sigma is 0.00149 RAW pds, i.e. "
          f"{0.00149/PDS_SPAN/6:.6f} overall.")
    json.dump(rows, open(WORK / "pdsham196.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
