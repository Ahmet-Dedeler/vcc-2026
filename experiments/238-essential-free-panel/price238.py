"""238 — the six-member composite of the October sources on every panel, on 228's exact routes.

`experiments/236-october-stack/price236.py`, unchanged in every constant and every route, so
the essential-free number is directly comparable to 236's published
**+0.0066…+0.0104 (pess) / +0.0057…+0.0081 (real)** and to the two arms AGENTS.md §8.1d cites
(228's magnitude arm −0.00702; 229's `ARM6` −0.00219).  The only differences are the input file
and the extra panels.

  `pds`   d_pds / 0.452 / 6 * 1.011                              (162; LOO sigma 0.001494)
  `mse`   197's inversion: d_raw_mse = -2.074495 * d_cos_board, then / 0.950162566 / 6,
          with `d_cos_board = OURS_COS * rel_cos * k_cos`, `k_cos` carried as a RANGE
          (1.000 on `h1union`, 1.384 on `cd4f034`) because 219 showed it is not identified
  `nmae`  188's identity on the `wald89` pool, / 0.5963 / 6, sign flipped (lower is better)
  `fid`   the shipped call set is prescribed, so the yield term is pinned at 1; / 0.30 / 6
  `jac`   only the 400-gene head can move; rescaled onto the board's raw `jac` 0.04096
  `reach` rule 9: its normalisation carries `n_real`, unreadable on any harness we hold, so it
          is quoted SEPARATELY and the composite is given WITH and WITHOUT it

Every number is reported RAW and NULL-CORRECTED (arm minus its own whole-row derangement at the
same dose), which is 233 §3.2's methodology and the correction 236 §3.1 showed is load-bearing.

🔴 The constants were calibrated on TODAY's panel and today's board.  These arms are measured on
synthetic panels that do not exist yet, so a composite here is a SIZE, not a board prediction.

    uv run python price238.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panels238 as P  # noqa: E402

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
CLEAN2 = ("orion_hct116", "orion_hek293t")
SPAN_PDS, CALIB, SIGMA_PDS = 0.452, 1.011, 0.001494
SPAN_MSE, DMSE_DCOS = 0.950162566, 2.074495
SPAN_NMAE = 0.5963
SPAN_FID, SPAN_JAC, SPAN_REACH = 0.30, 0.365, 0.895
OURS_COS = 0.09946902740184964
OURS_JAC_RAW, OURS_REACH_RAW = 0.04095656108311888, 0.27654544083789157
K_COS = (1.000, 1.384)
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089
ORDERING_NOISE = 0.004502          # §8.1b: sd on overall between digest-identical arms


def board_units(rel, dpds, dnm, dfid, rjac, rrch) -> dict:
    o_pds = dpds / SPAN_PDS / 6 * CALIB
    o_mse = {k: (DMSE_DCOS * OURS_COS * rel / 100 * k) / SPAN_MSE / 6 for k in K_COS}
    o_nmae = -dnm / SPAN_NMAE / 6
    o_fid = dfid / SPAN_FID / 6
    o_jac = (rjac * OURS_JAC_RAW) / SPAN_JAC / 6
    o_reach = (rrch * OURS_REACH_RAW) / SPAN_REACH / 6
    comp = {k: o_pds + v + o_nmae + o_fid + o_jac for k, v in o_mse.items()}
    return dict(o_pds=o_pds, o_mse=o_mse, o_nmae=o_nmae, o_fid=o_fid, o_jac=o_jac,
                o_reach=o_reach, composite_no_reach=comp,
                composite_with_reach={k: v + o_reach for k, v in comp.items()})


def main() -> None:
    rows = [r for r in json.loads((P.WORK / "six238.json").read_text())
            if not r.get("vacuous")]
    by = {(r["panel"], r["fold"], r["arm"]): r for r in rows}
    panels = [p for p in ("abc", "real", "pess", "efr", "efp")
              if any(r["panel"] == p for r in rows)]
    out: dict = {}

    for panel in panels:
        arms = sorted({r["arm"] for r in rows
                       if r["panel"] == panel and r["arm"] not in ("BASE", "null")})
        if not arms:
            continue
        print(f"\n{'=' * 108}\nPANEL {panel}\n{'=' * 108}")
        print(f"{'arm':<44}{'folds':>7}{'carry':>7}{'relcos%':>9}{'d_pds':>10}"
              f"{'d_nmae':>10}{'d_fid':>9}{'d_jac%':>9}{'d_reach%':>10}{'e=dM/dU':>10}")
        print("-" * 125)
        for arm in arms:
            for foldset, lab in ((CLEAN2, "CLEAN2"), (CANCER, "ALL3")):
                pairs = [(by[(panel, f, "BASE")], by[(panel, f, arm)]) for f in foldset
                         if (panel, f, arm) in by and (panel, f, "BASE") in by]
                if len(pairs) != len(foldset):
                    continue
                m = lambda fn: float(np.mean([fn(b, a) for b, a in pairs]))  # noqa: E731
                rel = m(lambda b, a: 100 * (a["cos_agg"] - b["cos_agg"]) / abs(b["cos_agg"]))
                dpds = m(lambda b, a: a["pds"] - b["pds"])
                dnm = m(lambda b, a: a["nmae_pools"]["wald89"] - b["nmae_pools"]["wald89"])
                dfid = m(lambda b, a: a["fid"] - b["fid"])
                rjac = m(lambda b, a: a["jac"]["89"] / b["jac"]["89"] - 1.0)
                rrch = m(lambda b, a: a["reach"] / b["reach"] - 1.0)
                du = m(lambda b, a: a["u_over_d"] - b["u_over_d"])
                dm = m(lambda b, a: a["m_over_d"] - b["m_over_d"])
                e = dm / du if du != 0 else float("nan")
                bu = board_units(rel, dpds, dnm, dfid, rjac, rrch)
                pos = sum(1 for b, a in pairs if a["pds"] > b["pds"])
                carry = pairs[0][1].get("carrying_rows")
                out[f"{panel}|{arm}|{lab}"] = dict(
                    rel_cos=rel, d_pds=dpds, d_nmae=dnm, d_fid=dfid, rel_jac=rjac,
                    rel_reach=rrch, dU=du, dM=dm, e=e, carrying_rows=carry,
                    carrying_essential=pairs[0][1].get("carrying_targets_essential"),
                    pds_positive=f"{pos}/{len(pairs)}", **bu)
                if lab == "CLEAN2":
                    c = bu["composite_no_reach"]
                    cr = bu["composite_with_reach"]
                    print(f"{arm:<44}{lab:>7}{str(carry):>7}{rel:+9.2f}{dpds:+10.5f}"
                          f"{dnm:+10.5f}{dfid:+9.5f}{100 * rjac:+9.2f}"
                          f"{100 * rrch:+10.2f}{e:+10.3f}")
                    print(f"{'    -> board units':<44}{'':>7}"
                          f"  pds {bu['o_pds']:+.5f}  mse {bu['o_mse'][1.0]:+.5f}.."
                          f"{bu['o_mse'][1.384]:+.5f}  nmae {bu['o_nmae']:+.5f}  "
                          f"fid {bu['o_fid']:+.5f}  jac {bu['o_jac']:+.5f}  "
                          f"[reach {bu['o_reach']:+.5f}]")
                    print(f"{'    -> COMPOSITE':<44}{'':>7}"
                          f"  {c[1.0]:+.5f} .. {c[1.384]:+.5f}   (with reach "
                          f"{cr[1.0]:+.5f} .. {cr[1.384]:+.5f})   pds+ {pos}/{len(pairs)}"
                          f"   gates: shape {GATE_SHAPE:+.4f} / amp {GATE_AMP:+.4f}")

    # ---------------------------------------------------------------- null-corrected price
    print(f"\n{'=' * 108}\nNULL-CORRECTED (arm minus its OWN derangement at the same dose)"
          f"\n{'=' * 108}")
    print(f"{'panel  arm':<52}{'relcos%':>9}{'d_pds':>10}{'d_nmae':>10}{'d_fid':>9}"
          f"{'e':>9}{'COMPOSITE (no reach)':>28}")
    for key, v in list(out.items()):
        panel, arm, lab = key.split("|")
        if lab != "CLEAN2" or "_shuf@" in arm:
            continue
        base, _, dose = arm.rpartition("@")
        s = out.get(f"{panel}|{base}_shuf@{dose}|CLEAN2")
        if s is None:
            continue
        rel = v["rel_cos"] - s["rel_cos"]
        dpds = v["d_pds"] - s["d_pds"]
        dnm = v["d_nmae"] - s["d_nmae"]
        dfid = v["d_fid"] - s["d_fid"]
        rjac = v["rel_jac"] - s["rel_jac"]
        rrch = v["rel_reach"] - s["rel_reach"]
        du, dm = v["dU"] - s["dU"], v["dM"] - s["dM"]
        bu = board_units(rel, dpds, dnm, dfid, rjac, rrch)
        c = bu["composite_no_reach"]
        out[f"{key}|NULLCORR"] = dict(rel_cos=rel, d_pds=dpds, d_nmae=dnm, d_fid=dfid,
                                      rel_jac=rjac, rel_reach=rrch, dU=du, dM=dm,
                                      e=(dm / du if du else float("nan")),
                                      carrying_rows=v["carrying_rows"],
                                      carrying_essential=v["carrying_essential"], **bu)
        print(f"{panel + '  ' + arm:<52}{rel:+9.2f}{dpds:+10.5f}{dnm:+10.5f}{dfid:+9.5f}"
              f"{(dm / du if du else float('nan')):+9.3f}"
              f"{c[1.0]:+15.5f}..{c[1.384]:+.5f}")

    print("\nLEAVE-TOP-3-OUT RETENTION (falsifier)")
    for r in rows:
        if "leave_top3" in r:
            lt = r["leave_top3"]
            print(f"  {r['panel']:5s} {r['fold']:14s} {r['arm']:<46} "
                  f"cos {lt['cos_gain_full']:+.6f} -> {lt['cos_gain_minus_top3']:+.6f} "
                  f"({100 * lt['cos_retained']:6.1f}%)   "
                  f"pds {lt['pds_gain_full']:+.6f} -> {lt['pds_gain_minus_top3']:+.6f}")

    (P.WORK / "price238.json").write_text(json.dumps(out, indent=2, default=float))
    print(f"\n-> {P.WORK / 'price238.json'}")


if __name__ == "__main__":
    main()
