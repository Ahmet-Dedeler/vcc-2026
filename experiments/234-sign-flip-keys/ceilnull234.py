"""234 — two control numbers the headline needs.

(1) THE CEILING NULL.  Phase 2's in-fold calibrated ceiling flips every bin where
    `P(correct | key) < 0.5`, with the bins fitted on the fold being scored.  A PURE-NOISE
    key fitted the same way also produces a positive number, purely from binning noise.
    This measures that floor at 20 / 50 / 200 bins, so the real keys' ceilings can be read
    against it instead of against zero.

(2) THE EFFICIENCY OF A REAL KEY AGAINST 232's ORACLE.  232's exchange rate was fitted on an
    oracle that flips ONLY wrong signs (every flip is a gain).  A key achieves a given net
    `Dq_w` by making both right and wrong flips, so it moves far more magnitude for the same
    net gain.  Comparing 234's `combo_pool` rows to 232's `oracle` rows at matched `Dq_w`
    gives the conversion factor between "oracle `Dq_w`" and "key `Dq_w`", i.e. how much
    higher the real bar is.

    uv run python ceilnull234.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from common234 import CANCER, GATE, QW_NEEDED, RATE_PER_QW, SEED, WORK  # noqa: E402
from phase0_234 import six_member  # noqa: E402
from screen234 import Fold, calibrated_ceiling  # noqa: E402

NBINS = (20, 50, 200)
NDRAW = 30


def main() -> None:
    rng = np.random.default_rng(SEED)
    folds = [Fold(h) for h in CANCER]
    out = {}

    print("\n(1) CEILING NULL — a PURE-NOISE key, calibrated in-fold exactly like a real one")
    print(f"{'bins':>6}{'null mean':>13}{'null p95':>12}{'null max':>12}"
          f"{'best real key':>16}{'combo key':>12}")
    print("-" * 71)
    screen = json.loads((WORK / "screen234.json").read_text())
    honest = json.loads((WORK / "honest234.json").read_text())
    for nb in NBINS:
        vals = []
        for d in range(NDRAW):
            v = [calibrated_ceiling(f, rng.random(f.truth.shape).astype(np.float32), nb)[0]
                 for f in folds]
            vals.append(float(np.mean(v)))
        v = np.asarray(vals)
        best_real = max(r[f"ceiling_infold_{nb}"] for r in screen["phase2"])
        combo = honest["ceiling_combo"][f"infold_{nb}"]
        out[f"null_{nb}"] = dict(mean=float(v.mean()), p95=float(np.percentile(v, 95)),
                                 max=float(v.max()), sd=float(v.std(ddof=1)),
                                 best_real_key=best_real, combo_key=combo)
        print(f"{nb:>6d}{v.mean():>+13.5f}{np.percentile(v,95):>+12.5f}{v.max():>+12.5f}"
              f"{best_real:>+16.5f}{combo:>+12.5f}")

    print("\n(2) THE EFFICIENCY OF A REAL KEY vs 232's TRUTH-SIGN ORACLE, at matched d q_w")
    rows232 = json.loads(
        (ROOT / "experiments/232-sign-vote/work/arm232.json").read_text())
    R2 = {(r["fold"], r["arm"]): r for r in rows232}
    base2 = [R2[(f, "BASE")] for f in CANCER]
    orc = {}
    for arm in sorted({r["arm"] for r in rows232 if "|oracle@" in r["arm"]}):
        if not all((f, arm) in R2 for f in CANCER):
            continue
        p = six_member(base2, [R2[(f, arm)] for f in CANCER], 1.384)
        orc[arm] = p
    rows234 = json.loads((WORK / "arm234.json").read_text())
    R4 = {(r["fold"], r["arm"]): r for r in rows234}
    base4 = [R4[(f, "BASE")] for f in CANCER]

    print(f"{'arm':<26}{'d q_w':>10}{'composite':>12}{'per unit q_w':>15}{'vs oracle':>12}")
    print("-" * 75)
    ref = RATE_PER_QW["k1.384"]

    # (2a) 232's rate is a TWO-SIDED fit.  Split it: the response is asymmetric, not noisy.
    def sl(ps):
        x = np.array([a for a, _ in ps])
        y = np.array([b for _, b in ps])
        return (float((x * y).sum() / (x * x).sum()),
                float(np.corrcoef(x, y)[0, 1]) if len(x) > 2 else float("nan"))
    allp = []
    for arm in sorted({r["arm"] for r in rows232
                       if "oracle@" in r["arm"]}):                # includes ANTIoracle
        if not all((f, arm) in R2 for f in CANCER):
            continue
        p = six_member(base2, [R2[(f, arm)] for f in CANCER], 1.384)
        allp.append((p["d_qw"], p["composite"]))
    posp = [p for p in allp if p[0] > 0]
    negp = [p for p in allp if p[0] < 0]
    sides = {}
    for nm, ps in (("both_14", allp), ("oracle_positive_7", posp),
                   ("antioracle_negative_7", negp)):
        s, r = sl(ps)
        sides[nm] = dict(slope=s, per_0p01=s * 0.01, r=r, n=len(ps),
                         qw_needed_for_gate=GATE / s)
    out["exchange_rate_sides"] = sides
    print("  232's exchange rate, SPLIT BY SIDE (its published +0.003990 pools both):")
    for nm, v in sides.items():
        print(f"    {nm:<24} {v['per_0p01']:+.6f} per +0.01 q_w   r {v['r']:+.4f}   "
              f"gate needs d q_w +{v['qw_needed_for_gate']:.5f}")
    print("  -> the response is ASYMMETRIC (r = 1.000 on each side separately): losing sign "
          "accuracy costs\n     1.88x what gaining it pays.  For pricing a POSITIVE sign arm "
          "the one-sided rate is the right one.\n")
    for arm in sorted(orc, key=lambda a: orc[a]["d_qw"]):
        p = orc[arm]
        if p["d_qw"] <= 0:
            continue
        print(f"{'232 ' + arm.split('|')[-1]:<26}{p['d_qw']:>+10.5f}"
              f"{p['composite']:>+12.5f}{p['composite']/p['d_qw']:>15.4f}"
              f"{p['composite']/p['d_qw']/ref:>12.2f}x")
    key_pts = []
    for arm in sorted({r["arm"] for r in rows234 if r["arm"].startswith("combo_pool@")}):
        if not all((f, arm) in R4 for f in CANCER):
            continue
        p = six_member(base4, [R4[(f, arm)] for f in CANCER], 1.384)
        if p["d_qw"] <= 1e-5:
            continue
        key_pts.append((p["d_qw"], p["composite"]))
        print(f"{'234 ' + arm:<26}{p['d_qw']:>+10.5f}{p['composite']:>+12.5f}"
              f"{p['composite']/p['d_qw']:>15.4f}{p['composite']/p['d_qw']/ref:>12.2f}x")

    def slope(pts):
        x = np.array([a for a, _ in pts])
        y = np.array([b for _, b in pts])
        return float((x * y).sum() / (x * x).sum())

    hi = max(a for a, _ in key_pts)
    orc_pts = [(p["d_qw"], p["composite"]) for p in orc.values() if 0 < p["d_qw"]]
    orc_matched = [(a, b) for a, b in orc_pts if a <= 3.0 * hi]
    s_key, s_orc_m, s_orc_all = slope(key_pts), slope(orc_matched), slope(orc_pts)
    out["key_vs_oracle_efficiency"] = dict(
        slope_key=s_key, slope_oracle_matched_range=s_orc_m, slope_oracle_all=s_orc_all,
        published_rate=ref, eff_vs_matched=s_key / s_orc_m, eff_vs_published=s_key / ref,
        effective_qw_bar_matched=QW_NEEDED["k1.384"] * s_orc_m / s_key,
        effective_qw_bar_published=QW_NEEDED["k1.384"] * ref / s_key,
        key_points=key_pts, oracle_points=orc_pts, matched_range_hi=3.0 * hi)
    print(f"\n  through-origin slope, composite per unit d q_w:")
    print(f"    232 oracle, ALL 7 positive doses            {s_orc_all:.4f}")
    print(f"    232 oracle, matched d q_w range (<= {3*hi:.4f})  {s_orc_m:.4f}")
    print(f"    234 key (combo_pool, pool-confined)         {s_key:.4f}")
    print(f"  -> a real key realises {s_key/s_orc_m:.2f}x of the oracle at MATCHED dose, "
          f"{s_key/ref:.2f}x of 232's published rate ({ref:.4f}).")
    print(f"  -> the EFFECTIVE bar on a key's d q_w is "
          f"+{QW_NEEDED['k1.384']*s_orc_m/s_key:.5f} (matched) .. "
          f"+{QW_NEEDED['k1.384']*ref/s_key:.5f} (published), not +{QW_NEEDED['k1.384']:.5f}.")
    (WORK / "ceilnull234.json").write_text(json.dumps(out, indent=1))
    print(f"\n-> {WORK/'ceilnull234.json'}")


if __name__ == "__main__":
    main()
