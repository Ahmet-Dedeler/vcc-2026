"""236 Phase 2 — the six-member composite of the October stack, on 228's EXACT routes.

Every constant and every route below is copied from `experiments/228-theta-gain-plolo/
price228.py`, unchanged, so this number is comparable to the two the mandate cites: 228's
magnitude arm (**-0.00702**) and 229's `ARM6` (**-0.00219** realised against a +0.00268 rate
card).  The only thing that differs is the arm and the panel.

  `pds`   d_pds / 0.452 / 6 * 1.011                              (162; LOO sigma 0.001494)
  `mse`   197's inversion: d_raw_mse = -2.074495 * d_cos_board, then / 0.950162566 / 6,
          with `d_cos_board = OURS_COS * rel_cos * k_cos` and `k_cos` carried as a RANGE
          (1.000 on `h1union`, 1.384 on `cd4f034`) because it is not identified.
  `nmae`  188's exact identity on the `wald89` pool, / 0.5963 / 6, sign flipped (lower better)
  `fid`   the shipped call set is prescribed (400 + 2,200), so the yield term is pinned at 1
          and `fid` is the call set's sign precision alone; / 0.30 / 6
  `jac`   only the 400-gene head can move; relative delta rescaled onto the board's own raw
          `jac` 0.04096; / 0.365 / 6
  `reach` rule 9: its normalisation carries `n_real`, which is unreadable on any harness we
          hold, so it is quoted SEPARATELY and the composite is given WITH and WITHOUT it.

🔴 READ THIS BEFORE QUOTING THE COMPOSITE.  These constants were calibrated on TODAY's panel
and today's board.  The arm is measured on a DISJOINT synthetic panel that does not exist yet.
So the composite is a SIZE, not a board prediction, and 233 §6 is the reason: Arc SELECTS its
panel and 199's are random draws, which makes every number here a lower bound of unknown
slack.  Nothing here is submittable before 2026-10-22.

    uv run python price236.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

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


def main() -> None:
    rows = json.loads((WORK / "six236.json").read_text())
    by = {(r["panel"], r["fold"], r["arm"]): r for r in rows}
    panels = sorted({r["panel"] for r in rows})
    out: dict = {}

    for panel in panels:
        arms = sorted({r["arm"] for r in rows if r["panel"] == panel and r["arm"] != "BASE"})
        print(f"\n{'=' * 100}\nPANEL {panel}\n{'=' * 100}")
        print(f"{'arm':<44}{'folds':>6}{'relcos%':>9}{'d_pds':>10}{'d_nmae':>10}"
              f"{'d_fid':>9}{'d_jac%':>9}{'d_reach%':>10}{'e=dM/dU':>10}")
        print("-" * 107)
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

                o_pds = dpds / SPAN_PDS / 6 * CALIB
                o_mse = {k: (DMSE_DCOS * OURS_COS * rel / 100 * k) / SPAN_MSE / 6
                         for k in K_COS}
                o_nmae = -dnm / SPAN_NMAE / 6
                o_fid = dfid / SPAN_FID / 6
                o_jac = (rjac * OURS_JAC_RAW) / SPAN_JAC / 6
                o_reach = (rrch * OURS_REACH_RAW) / SPAN_REACH / 6
                comp = {k: o_pds + v + o_nmae + o_fid + o_jac for k, v in o_mse.items()}
                comp_r = {k: v + o_reach for k, v in comp.items()}
                pos = sum(1 for b, a in pairs if a["pds"] > b["pds"])

                out[f"{panel}|{arm}|{lab}"] = dict(
                    rel_cos=rel, d_pds=dpds, d_nmae=dnm, d_fid=dfid, rel_jac=rjac,
                    rel_reach=rrch, dU=du, dM=dm, e=e, pds_positive=f"{pos}/{len(pairs)}",
                    o_pds=o_pds, o_mse=o_mse, o_nmae=o_nmae, o_fid=o_fid, o_jac=o_jac,
                    o_reach=o_reach, composite_no_reach=comp, composite_with_reach=comp_r)
                if lab == "CLEAN2":
                    print(f"{arm:<44}{lab:>6}{rel:+9.2f}{dpds:+10.5f}{dnm:+10.5f}"
                          f"{dfid:+9.5f}{100*rjac:+9.2f}{100*rrch:+10.2f}{e:+10.3f}")
                    print(f"{'    -> board units':<44}{'':>6}"
                          f"  pds {o_pds:+.5f}  mse {o_mse[1.0]:+.5f}..{o_mse[1.384]:+.5f}"
                          f"  nmae {o_nmae:+.5f}  fid {o_fid:+.5f}  jac {o_jac:+.5f}"
                          f"  [reach {o_reach:+.5f}]")
                    print(f"{'    -> COMPOSITE':<44}{'':>6}"
                          f"  {comp[1.0]:+.5f} .. {comp[1.384]:+.5f}   (with reach "
                          f"{comp_r[1.0]:+.5f} .. {comp_r[1.384]:+.5f})   pds+ {pos}/{len(pairs)}"
                          f"   gates: shape {GATE_SHAPE:+.4f} / amplitude {GATE_AMP:+.4f}")

    # ---------------------------------------------------------------- null-corrected price
    # 233 §3.2's methodology: every number is CREDIT ABOVE A MEASURED NULL, not a raw delta.
    # The null here is the arm's OWN whole-row derangement at the SAME dose -- it holds
    # magnitude fixed, so anything the derangement also delivers is a dose artefact and not
    # the source.  This matters: on the realistic panel the derangement reproduces the whole
    # `nmae` credit, which is why the raw d_nmae disagrees between panels and the corrected
    # one does not.
    print(f"\n{'=' * 100}\nNULL-CORRECTED (arm minus its OWN derangement at the same dose)\n{'=' * 100}")
    print(f"{'panel  arm':<52}{'relcos%':>9}{'d_pds':>10}{'d_nmae':>10}{'d_fid':>9}"
          f"{'e':>9}{'COMPOSITE (no reach)':>26}")
    for key, v in list(out.items()):
        panel, arm, lab = key.split("|")
        if lab != "CLEAN2" or arm.endswith("_shuf@0.34"):
            continue
        base, _, dose = arm.rpartition("@")
        sk = f"{panel}|{base}_shuf@{dose}|CLEAN2"
        s = out.get(sk)
        if s is None:
            continue
        rel = v["rel_cos"] - s["rel_cos"]
        dpds = v["d_pds"] - s["d_pds"]
        dnm = v["d_nmae"] - s["d_nmae"]
        dfid = v["d_fid"] - s["d_fid"]
        rjac = v["rel_jac"] - s["rel_jac"]
        rrch = v["rel_reach"] - s["rel_reach"]
        du, dm = v["dU"] - s["dU"], v["dM"] - s["dM"]
        o_pds = dpds / SPAN_PDS / 6 * CALIB
        o_mse = {k: (DMSE_DCOS * OURS_COS * rel / 100 * k) / SPAN_MSE / 6 for k in K_COS}
        o_nmae = -dnm / SPAN_NMAE / 6
        o_fid = dfid / SPAN_FID / 6
        o_jac = (rjac * OURS_JAC_RAW) / SPAN_JAC / 6
        o_reach = (rrch * OURS_REACH_RAW) / SPAN_REACH / 6
        comp = {k: o_pds + val + o_nmae + o_fid + o_jac for k, val in o_mse.items()}
        out[f"{key}|NULLCORR"] = dict(
            rel_cos=rel, d_pds=dpds, d_nmae=dnm, d_fid=dfid, rel_jac=rjac, rel_reach=rrch,
            dU=du, dM=dm, e=(dm / du if du else float("nan")), o_pds=o_pds, o_mse=o_mse,
            o_nmae=o_nmae, o_fid=o_fid, o_jac=o_jac, o_reach=o_reach,
            composite_no_reach=comp,
            composite_with_reach={k: val + o_reach for k, val in comp.items()})
        print(f"{panel + '  ' + arm:<52}{rel:+9.2f}{dpds:+10.5f}{dnm:+10.5f}{dfid:+9.5f}"
              f"{(dm / du if du else float('nan')):+9.3f}"
              f"{comp[1.0]:+13.5f}..{comp[1.384]:+.5f}")

    print("\nLEAVE-TOP-3-OUT RETENTION (falsifier 4)")
    for r in rows:
        if "leave_top3" in r:
            lt = r["leave_top3"]
            print(f"  {r['panel']:5s} {r['fold']:14s} {r['arm']:<46} "
                  f"cos {lt['cos_gain_full']:+.6f} -> {lt['cos_gain_minus_top3']:+.6f} "
                  f"({100*lt['cos_retained']:6.1f}%)   "
                  f"pds {lt['pds_gain_full']:+.6f} -> {lt['pds_gain_minus_top3']:+.6f} "
                  f"({100*lt['pds_retained']:7.1f}%)")

    (WORK / "price236.json").write_text(json.dumps(out, indent=2, default=float))
    print(f"\n-> {WORK / 'price236.json'}")


if __name__ == "__main__":
    main()
