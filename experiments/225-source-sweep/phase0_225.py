"""225 Phase 0 — MANDATORY known-answer rows.  Nothing after this counts if it fails.

Reproduces, through 225's own harness object (`harness225.Fold`), rather than by re-reading
223's JSON:

  KA-1  the three baselines `f0` / `cd4@0.34` / `base222` == 218's `resid218` +
        `compose218` rows, on every fold, to 1e-9.  (223's KA-1/2/3.)
  KA-2  223's HEADLINE: `dld1_lowmoi` @ g=0.35 marginal over base222, 3 cancer folds
        -> rel Δmacro-cos +1.405%, Δpds +0.00256, six-member price +0.00228.
  KA-3  two of 223's seven deaths, reproduced with their own numbers:
          `telohaec`   @ g=0.28  -> rel Δcos -0.410%, Δpds -0.00145, price -0.00120
          `cmap_xpr`   @ g=0.28  -> rel Δcos +1.048%, Δpds -0.00138, price -0.00077
        (a third, `kolf2` @ 0.17 -> -0.00174, is run as a bonus row)
  KA-4  the silent-zero detector still works: `dld1` under the OLD median-over-all-rows
        rule must be EXACTLY bit-identical to base222, and `kolf2` (91%) must NOT be.

    uv run python phase0_225.py
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from harness225 import (CANCER, GATE_SHAPE_ONLY, WORK, Fold, price,  # noqa: E402
                        solve_fraction, source_residual)

TOL_COS_PCT = 0.05      # 223's own instrument resolution is +-0.20% relative (204)
TOL_PDS = 5e-5
TOL_PRICE = 5e-5

# 223 RESULTS.md tables 4.1 / 4.2, three cancer folds, marginal over base222
EXPECT = {
    ("dld1_lowmoi", 0.35): dict(rel_cos=+1.405, d_pds=+0.00256, price=+0.00228),
    ("dld1_lowmoi", 0.28): dict(rel_cos=+1.403, d_pds=+0.00250, price=+0.00224),
    ("telohaec", 0.28): dict(rel_cos=-0.410, d_pds=-0.00145, price=-0.00120),
    ("cmap_xpr", 0.28): dict(rel_cos=+1.048, d_pds=-0.00138, price=-0.00077),
    ("kolf2", 0.17): dict(rel_cos=+0.544, d_pds=-0.00256, price=-0.00174),
}


def main() -> None:
    ref_resid = {(r["fold"], r["arm"]): r for r in json.loads(
        (Path(__file__).resolve().parents[1]
         / "218-independent-donors/work/resid218.json").read_text())}
    ref_comp = {(r["fold"], r["arm"]): r for r in json.loads(
        (Path(__file__).resolve().parents[1]
         / "218-independent-donors/work/compose218.json").read_text())}

    ka: list[dict] = []
    arms: list[dict] = []
    for held in CANCER:
        t0 = time.time()
        F = Fold(held)
        print(f"[{held}] fold built in {time.time() - t0:.0f}s  "
              f"frac_cd4 {F.frac_cd4:.5f} frac_union {F.frac_union:.5f}", flush=True)

        # -------------------------------------------------------------- KA-1 baselines
        for arm, refmap, refarm in (("f0", ref_resid, "f0"),
                                    ("cd4", ref_resid, "cd4_gn0.34"),
                                    ("base222", ref_comp, "cd40.34+h1PAIRCARRY0.28")):
            r = refmap[(held, refarm)]
            dc = abs(r["cos_agg"] - F.scored[arm]["cos_agg"])
            dp = abs(r["pds"] - F.scored[arm]["pds"])
            ok = dc < 1e-9 and dp < 1e-9
            ka.append(dict(gate=f"KA-1 {held}/{arm}", ref=f"218 {refarm}",
                           d_cos_agg=dc, d_pds=dp, status="PASS" if ok else "FAIL"))
            assert ok, f"KA-1 FAILED {held}/{arm}: dcos {dc:.3e} dpds {dp:.3e}"
            print(f"    KA-1 PASS {arm} == 218 {refarm} "
                  f"(dcos {dc:.1e} dpds {dp:.1e})", flush=True)

        # ------------------------------------------------------ KA-2/3 the five re-prices
        for (sname, g) in EXPECT:
            if sname == held:
                continue
            S = source_residual(sname, F.targets, "gene_norm")
            a = F.arm(S, g)
            a.update(fold=held, source=sname)
            arms.append(a)
            print(f"    [{held}/{sname}@{g}] relcos {a['rel_cos_pct']:+.3f}% "
                  f"dpds {a['d_pds']:+.5f} carry {a['audit']['carrying_rows']}"
                  f"/{a['audit']['n_rows']}", flush=True)

        # ------------------------------------------------------------- KA-4 silent zero
        for sname, expect_void in (("dld1_lowmoi", True), ("kolf2", False)):
            if sname == held:
                continue
            S = source_residual(sname, F.targets, "gene_norm")
            fr, _ = solve_fraction(F.whitened, S, 0.28)     # the OLD production rule
            void = fr == 0.0
            ok = void == expect_void
            ka.append(dict(gate=f"KA-4 {held}/{sname} OLD-rule void?",
                           expected=expect_void, measured=void, fraction=float(fr),
                           status="PASS" if ok else "FAIL"))
            assert ok, f"KA-4 FAILED {held}/{sname}: void={void} expected {expect_void}"
        print(f"    KA-4 PASS silent-zero detector live on {held}", flush=True)
        del F

    # ------------------------------------------------------- aggregate the cancer folds
    agg = {}
    for (sname, g), exp in EXPECT.items():
        sel = [a for a in arms if a["source"] == sname and a["g"] == g]
        if not sel:
            continue
        rel = float(np.mean([a["rel_cos_pct"] for a in sel]))
        dcos = float(np.mean([a["d_cos"] for a in sel]))
        dpds = float(np.mean([a["d_pds"] for a in sel]))
        p = price(dcos, dpds)
        ok = (abs(rel - exp["rel_cos"]) < TOL_COS_PCT
              and abs(dpds - exp["d_pds"]) < TOL_PDS
              and abs(p["total"] - exp["price"]) < TOL_PRICE)
        agg[f"{sname}@{g}"] = dict(n_folds=len(sel), rel_cos_pct=rel, d_cos=dcos,
                                   d_pds=dpds, price_k1=p["total"],
                                   expected=exp, status="PASS" if ok else "FAIL")
        ka.append(dict(gate=f"KA-{'2' if sname == 'dld1_lowmoi' else '3'} "
                            f"223 {sname}@{g}", **agg[f"{sname}@{g}"]))
        print(f"  {sname}@{g:<5g} rel {rel:+.3f}% (exp {exp['rel_cos']:+.3f}) "
              f"dpds {dpds:+.5f} (exp {exp['d_pds']:+.5f}) "
              f"price {p['total']:+.5f} (exp {exp['price']:+.5f})  "
              f"{'PASS' if ok else 'FAIL'}", flush=True)

    (WORK / "phase0_225.json").write_text(json.dumps(
        dict(ka=ka, arms=arms, agg=agg, gate_shape_only=GATE_SHAPE_ONLY), indent=2))
    fails = [k for k in ka if k["status"] != "PASS"]
    print(f"\n{len(ka)} gates, {len(fails)} FAIL -> {WORK / 'phase0_225.json'}")
    if fails:
        raise SystemExit("PHASE 0 FAILED — the harness is wrong, nothing after it counts:\n"
                         + json.dumps(fails, indent=2))


if __name__ == "__main__":
    main()
