"""207 — known-answer gates, the dose-response table, and the three maximisers.

PRIMARY ENDPOINT: `cos_agg` = sqrt(1 - sse) on the four PLOLO folds, the energy-weighted
aggregate cosine.  That is the same functional 197 inverts off the board
(`cos = X / (sqrt(SigmaD) * mean_c sqrt(A_c))`), evaluated on a surface where the truth exists.

KNOWN-ANSWER GATES (rule 10's corollary -- abort if they fail):
  KA-1  a14 -> a15 (the shipped CD4 arm) must reproduce 162's PLOLO prediction of
        +0.006656 SCALED `pds` over four folds.  This is the row that proves the dose axis
        here is the same axis production moved along.
  KA-2  a23 - a15 (the mass fix at the shipped dose) must reproduce 162's +0.004014 scaled.
  KA-3  the a14 -> a15 relative cosine gain must be near the board's +6.8% (197 section 1.2).
  KA-4  `f0.1` must equal `a15` exactly (the same configuration reached by two code paths).

    uv run python analyse207.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
SPAN_PDS, CALIB = 0.452, 1.011
BOARD_KA = dict(a15_pds_scaled=0.006656, a23_pds_scaled=0.004014, cos_rel=0.0683)
NB = 4000


def load(*names) -> dict:
    by: dict[tuple[str, str], dict] = {}
    for n in names:
        p = WORK / n
        if p.exists():
            for r in json.loads(p.read_text()):
                by[(r["fold"], r["arm"])] = r
    return by


def mean_over_folds(by, arm, key, folds=FOLDS):
    v = [by[(f, arm)][key] for f in folds if (f, arm) in by]
    return float(np.mean(v)) if len(v) == len(folds) else None


def paired_boot(by, arm, base, key, folds=FOLDS, seed=20260903):
    """Paired TARGET bootstrap of an arm-to-arm difference, 204's construction."""
    rng = np.random.default_rng(seed)
    diffs = []
    for f in folds:
        a = np.asarray(by[(f, arm)][key])
        b = np.asarray(by[(f, base)][key])
        idx = rng.integers(0, len(a), (NB, len(a)))
        diffs.append(np.array([(a[i] - b[i]).mean() for i in idx]))
    d = np.mean(diffs, axis=0)
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)), float((d > 0).mean())


def row(by, arm, base=None, folds=FOLDS):
    out = {"arm": arm}
    for k in ("pds", "cos_agg", "cos", "a15_sgn_w", "a15_sgn",
              "median_row_norm_final", "median_row_norm_wcols"):
        out[k] = mean_over_folds(by, arm, k, folds)
    out["a15_pds"] = mean_over_folds(by, arm, "a15_pds", folds)
    a = by[(folds[0], arm)]["audit"]
    out["fraction"] = a.get("fraction")
    out["eff_f"] = a.get("effective_f")
    if base is not None and out["pds"] is not None:
        for k in ("pds", "cos_agg", "cos", "a15_sgn_w"):
            out["d_" + k] = out[k] - mean_over_folds(by, base, k, folds)
        lo, hi, p = paired_boot(by, arm, base, "_cos_per_target_a23", folds)
        out["d_cos_ci"], out["d_cos_p"] = (lo, hi), p
        lo, hi, p = paired_boot(by, arm, base, "_pds_per_target_a23", folds)
        out["d_pds_ci"], out["d_pds_p"] = (lo, hi), p
    return out


def table(by, arms, base, folds=FOLDS, title=""):
    print(f"\n{title}")
    print(f"{'arm':<14}{'fraction':>9}{'eff f':>8}{'cos_agg':>10}{'Δcos_agg':>11}"
          f"{'Δcos 95% CI':>24}{'P>0':>6}{'pds':>10}{'Δpds':>10}{'P>0':>6}"
          f"{'q_w':>9}{'Δq_w':>10}{'rownorm':>9}")
    print("-" * 140)
    rows = []
    for a in arms:
        if (folds[0], a) not in by:
            continue
        r = row(by, a, base, folds)
        if r["pds"] is None:
            continue
        rows.append(r)
        ci = r.get("d_cos_ci", (float("nan"),) * 2)
        print(f"{a:<14}{r['fraction']:>9.5f}{(r['eff_f'] or 0):>8.4f}"
              f"{r['cos_agg']:>10.6f}{r.get('d_cos_agg', 0):>+11.6f}"
              f"  [{ci[0]:+.6f},{ci[1]:+.6f}]{r.get('d_cos_p', float('nan')):>6.3f}"
              f"{r['pds']:>10.6f}{r.get('d_pds', 0):>+10.6f}"
              f"{r.get('d_pds_p', float('nan')):>6.3f}"
              f"{r['a15_sgn_w']:>9.6f}{r.get('d_a15_sgn_w', 0):>+10.6f}"
              f"{r['median_row_norm_wcols']:>9.4f}")
    return rows


def main() -> None:
    by = load("sweep207.json")
    doses = [0.0, 0.01, 0.02, 0.035, 0.05, 0.07, 0.10, 0.14, 0.20, 0.28, 0.40, 0.55, 0.80]
    cd4_arms = [f"f{v:g}" for v in doses]

    # ------------------------------------------------------------------ known-answer gates
    print("=" * 100)
    print("KNOWN-ANSWER GATES")
    print("=" * 100)
    gates = []
    d_pds = mean_over_folds(by, "a15", "pds") - mean_over_folds(by, "a14", "pds")
    ka1 = d_pds / SPAN_PDS
    gates.append(("KA-1  a14->a15 PLOLO Δpds scaled (162: +0.006656)",
                  BOARD_KA["a15_pds_scaled"], ka1,
                  abs(ka1 - BOARD_KA["a15_pds_scaled"]) < 0.0015))
    d23 = mean_over_folds(by, "a15", "a23_pds") - mean_over_folds(by, "a15", "a15_pds")
    ka2 = d23 / SPAN_PDS
    gates.append(("KA-2  a15->a23 mass-fix PLOLO Δpds scaled (162: +0.004014)",
                  BOARD_KA["a23_pds_scaled"], ka2,
                  abs(ka2 - BOARD_KA["a23_pds_scaled"]) < 0.0015))
    c14 = mean_over_folds(by, "a14", "cos_agg")
    c15 = mean_over_folds(by, "a15", "cos_agg")
    ka3 = (c15 - c14) / c14
    gates.append(("KA-3  a14->a15 relative cos_agg gain (board 197 §1.2: +6.8%)",
                  BOARD_KA["cos_rel"], ka3, abs(ka3 - BOARD_KA["cos_rel"]) < 0.03))
    same = all(abs(by[(f, "a15")]["cos_agg"] - by[(f, "f0.1")]["cos_agg"]) < 1e-12
               for f in FOLDS)
    gates.append(("KA-4  `a15` == `f0.1` bit-identical (two code paths)", 1.0,
                  float(same), same))
    for name, exp, got, ok in gates:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<62} expected {exp:+.6f}  "
              f"measured {got:+.6f}")
    if not all(g[3] for g in gates):
        print("\n*** A KNOWN-ANSWER GATE FAILED — nothing below is readable. ***")

    # --------------------------------------------------------------------- the dose sweep
    print("\n" + "=" * 100)
    print("PHASE 1 — CD4 DOSE SWEEP, baseline = f0 (CD4 OFF).  Primary endpoint: cos_agg")
    print("=" * 100)
    cd4_rows = table(by, cd4_arms, "f0", title="CD4 residual (099 `f` frame)")
    print(f"\n  a14 (062's fraction 0.05, f=0.0528) and a15 (shipped f=0.10) for reference:")
    table(by, ["a14", "a15"], "f0", title="")

    shuf = [f"shuf{v:g}" for v in doses if v > 0]
    table(by, shuf, "f0", title="FALSIFIER 1 — target-shuffled CD4 (row permutation)")
    table(by, [f"self{v:g}" for v in (0.05, 0.10, 0.20, 0.40)], "f0",
          title="FALSIFIER — 099's self-residual (base mixture's own un-whitened effect)")
    table(by, [f"kolf{v:g}" for v in (0.05, 0.10, 0.20, 0.40)], "f0",
          folds=("k562_gwps", "orion_hct116", "orion_hek293t"),
          title="FALSIFIER 2 — KOLF2 (unmatched lineage, 273/300) as the residual "
                "[3 folds; kolf2 fold excluded as leakage]")

    # ------------------------------------------------------------------- the maximisers
    print("\n" + "=" * 100)
    print("THE THREE MAXIMISERS — the whole hypothesis is that they differ")
    print("=" * 100)
    for key, label in (("cos_agg", "aggregate cosine (197's route) [PRIMARY]"),
                       ("cos", "macro per-target cosine"),
                       ("pds", "PLOLO pds (062's objective)"),
                       ("a15_sgn_w", "|t|-weighted q_w (187 §5)")):
        vals = {r["arm"]: r[key] for r in cd4_rows}
        best = max(vals, key=vals.get)
        print(f"  {label:<42} argmax = {best:<8} at f = "
              f"{float(best[1:]):<6.3g} value {vals[best]:.6f}   "
              f"(shipped f0.1 = {vals['f0.1']:.6f}, "
              f"{100*(vals[best]/vals['f0.1']-1):+.2f}%)")

    # per-fold detail at a few doses
    print("\n" + "=" * 100)
    print("PER-FOLD cos_agg (falsifier 2's homogeneity test: is the gain lineage-specific?)")
    print("=" * 100)
    print(f"{'arm':<10}" + "".join(f"{f:>16}" for f in FOLDS) + f"{'  rel gain vs f0':>18}")
    for a in cd4_arms:
        if (FOLDS[0], a) not in by:
            continue
        v = [by[(f, a)]["cos_agg"] for f in FOLDS]
        b = [by[(f, "f0")]["cos_agg"] for f in FOLDS]
        rel = [100 * (x / y - 1) for x, y in zip(v, b)]
        print(f"{a:<10}" + "".join(f"{x:>16.6f}" for x in v)
              + "   " + " ".join(f"{r:+6.1f}%" for r in rel))

    # per-context (render frame, NOT truth lineage -- stated explicitly)
    print("\n" + "=" * 100)
    print("PER-RENDER-FRAME cos_agg at the shipped dose vs f0 (A/B/C are RENDER frames here;")
    print("the truth is the held DONOR line, so this is NOT a lineage measurement)")
    print("=" * 100)
    for ctx in "ABC":
        v = float(np.mean([by[(f, "f0.1")][f"a23_{ctx}_cos_agg"] for f in FOLDS]))
        b = float(np.mean([by[(f, "f0")][f"a23_{ctx}_cos_agg"] for f in FOLDS]))
        print(f"  frame {ctx}: f0 {b:.6f} -> f0.1 {v:.6f}   {100*(v/b-1):+.2f}%")

    json.dump({"gates": [(n, e, g, bool(o)) for n, e, g, o in gates],
               "cd4_rows": cd4_rows},
              open(WORK / "analyse207.json", "w"), indent=2)
    print(f"\nwrote {WORK/'analyse207.json'}")


if __name__ == "__main__":
    main()
