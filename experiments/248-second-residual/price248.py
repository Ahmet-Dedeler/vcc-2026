"""248 -- the rider priced on ALL SIX members under TWO instruments, and the mandate's gate.

Inputs, all measured earlier in this experiment or read from 243's artifacts on disk:
  work/repro248-summary.json   PLOLO through the builder: rel cos_agg / macro, d pds, d q_w,
                               the falsifiers (SHUF x2 seeds, CD4ROWS, rowsplit, leave-top-k)
  work/six248.json             216's exact `nmae` identity + the delivered-amplitude readout
  243 work/invert243.json      h1union's MEASURED board cosine and its (A_c, B_c, Kc_c, cap_c)
  243 work/rule243.json        the class factors (k_cos, k_pds) and the DE-four rates per step
  [work/amp248.json]           the package-level (A_c, B_c) when a package exists (post-build)

(i)  243 section 4.4's CLASS RULE with the arm's own delivered amplitude:
       d cos      = k_cos x rel_cos_agg(PLOLO, 3 cancer folds, per-fold mean) x cos_h1union(measured)
       mse leg    = 197's exact identity at the arm's own terms (amplitude r = the arm's delivered
                    amplitude relative to h1union: A -> r^2 A, B -> r^2 B) minus h1union's raw
       pds leg    = k_pds x d pds(PLOLO) / 0.452 / 6
       DE-four    = de4_per_dcos x d cos   (the sign-propagation rate; h1union's own 0.239)
     classes: new-source (k_cos 0.94, band 0.94-1.00; k_pds 0.70) = a first source into rows the
     union residual does not yet occupy, which is DLD-1's case (67/300, 233 rows empty of any
     union residual; CD4 occupies 276 but the union residual is what the class was defined on);
     re-dose (k_cos 1.36; k_pds 0.95-1.15) reported as sensitivity.
(ii) 216's EXACT `nmae` identity for the nmae member, reported SEPARATELY with its two-for-two
     mis-sign record; composite (ii) = (i)'s mse + pds + (DE-four with the nmae share replaced by
     the identity's number).

GATE (mandate C): (i) >= 1.5 x 0.0019 = 0.00285 AND d pds >= 0 on 3/3 cancer folds AND
leave-top-3-out retention >= 50% on pds.

    uv run python price248.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E243 = ROOT / "experiments/243-transfer-factor"
sys.path.insert(0, str(E243))
from lib243 import CTX, SUMD, board_row  # noqa: E402

SPAN_MSE, SPAN_PDS, SPAN_NMAE = 0.950162566, 0.452, 0.60
GATE = 0.0019
GATE_X = 1.5
CLASSES = {
    # label: (k_cos, k_pds, note)
    "new-source class, point (243 §4.4)": (0.94, 0.70, "h1union's own k_cos / k_pds"),
    "new-source class, a15-tag k": (1.00, 0.70, "the a15-tag reading of the same point"),
    "new-source class, HIGH band (k_pds 1.15)": (1.00, 1.15, "the band's upper edge"),
    "re-dose class (if DLD-1's rows count as occupied)": (1.36, 1.05, "CD4 f-steps; sensitivity only"),
}
K_DE4 = {"h1union point": 0.2390541399231271, "low (a14->a15)": 0.17342622793046183,
         "high (cd4f034)": 0.42278887033065365}
NMAE_SHARE_OF_DE4 = 0.4324      # the h1union board step: nmae +0.005881 of the four's +0.013601 scaled


def raw_mse(per, cos, r=1.0):
    """197's identity at amplitude r x the given terms (A -> r^2 A, B -> r^2 B, cap fixed)."""
    num = float(np.mean([r * r * per[c]["A"] + per[c]["Kc"] - min(per[c]["cap"], r * r * per[c]["B"]) for c in CTX]))
    msa = float(np.mean([r * np.sqrt(per[c]["A"]) for c in CTX]))
    X = cos * np.sqrt(SUMD) * msa
    return 1.0 + (num - 2.0 * X) / SUMD


def main() -> None:
    S = json.loads((WORK / "repro248-summary.json").read_text())
    SIX = json.loads((WORK / "six248.json").read_text())["arms"]
    INV = json.loads((E243 / "work/invert243.json").read_text())
    RULE = json.loads((E243 / "work/rule243.json").read_text())
    inc = INV["arms"]["h1union"]
    cos_inc = inc["inv"]["cos"]
    raw_inc = inc["board"]["raw"]["mse"]
    per_inc = inc["per"]
    chk = raw_mse(per_inc, cos_inc)
    amp = json.loads((WORK / "amp248.json").read_text()) if (WORK / "amp248.json").exists() else None
    b_inc = board_row("h1union")
    lines = []

    def P(t=""):
        print(t, flush=True)
        lines.append(t)

    P(f"KA: 197's identity at h1union's own terms reproduces its board raw mse {chk:.7f} vs {raw_inc:.7f} "
      f"(|d| {abs(chk - raw_inc):.1e}) {'PASS' if abs(chk - raw_inc) < 1e-9 else 'FAIL'}")
    assert abs(chk - raw_inc) < 1e-9
    P(f"h1union measured board cosine {cos_inc:.6f} (243, streamed off the package; the floor would say "
      f"{inc['floor']:.6f}); board overall {b_inc['overall']:.7f}")
    P(f"243's classes: {json.dumps({k: round(v, 4) for k, v in RULE['classes']['a23']['k'].items()})}, "
      f"k_pds per step { {s: round(RULE['rates'][s]['k_pds_3c'], 3) for s in RULE['rates']} }, "
      f"DE-four per d cos { {s: round(RULE['rates'][s]['de4_per_dcos'], 3) for s in RULE['rates']} }")

    results = {}
    for arm in ("CUM1", "CUM4"):
        key = f"{arm}@0.35"
        a = S["arms"][key]
        six_name = {"CUM1": "CUM1 dld1@0.35 <- RIDER", "CUM4": "CUM4 four@0.35"}[arm]
        six = SIX[six_name]
        conc = S["conc"][arm]["k"]["3"]
        rel_agg = a["rel_agg"]
        rel_macro = a["rel_macro"]
        d_pds = a["d_pds"]
        # delivered amplitude: the package's own (A_c, B_c) if built, else the effect-level read
        if amp is not None and arm in amp:
            r_amp = amp[arm]["amplitude_ratio"]
            amp_src = "PACKAGE (A_c streamed off the h5ad, 243's totals243)"
            per_arm = amp[arm]["per"]
        else:
            r_amp = six["delivered_norm_ratio"] / SIX["h1union g=0.28 <- INCUMBENT"]["delivered_norm_ratio"]
            amp_src = "EFFECT-LEVEL (216's delivered median row norm on frame A, no package yet)"
            per_arm = None
        u2_label = "dld1_lowmoi" if arm == "CUM1" else "225's four (dld1, jiang5, igvf, myllia)"
        P(f"\n{'=' * 110}\n{arm}: union 2 = {u2_label} @0.35 on top of the shipped "
          f"h1union (cd4@0.34 + H1 pair@0.28)\n{'=' * 110}")
        P(f"PLOLO, 3 cancer folds, through the builder: rel cos_agg {100*rel_agg:+.3f}%  rel macro {100*rel_macro:+.3f}%  "
          f"d pds {d_pds:+.6f} ({a['pds_pos']}/3 folds positive)  d q_w {a['d_qw']:+.5f}")
        P(f"  per fold: " + "; ".join(f"{f}: agg {100*v['rel_agg']:+.3f}% macro {100*v['rel_macro']:+.3f}% dpds {v['d_pds']:+.6f}"
                                       for f, v in a["per_fold"].items()))
        k2 = a.get("kolf2")
        if k2:
            P(f"  kolf2 (apart): rel macro {100*k2['rel_macro']:+.3f}%  rel agg {100*k2['rel_agg']:+.3f}%  d pds {k2['d_pds']:+.6f}")
        P(f"  falsifiers @0.35: SHUF20726 macro {100*S['arms'][f'{arm}SHUF20726@0.35']['rel_macro']:+.3f}% "
          f"(dpds {S['arms'][f'{arm}SHUF20726@0.35']['d_pds']:+.6f}); SHUF771103 macro "
          f"{100*S['arms'][f'{arm}SHUF771103@0.35']['rel_macro']:+.3f}% (dpds {S['arms'][f'{arm}SHUF771103@0.35']['d_pds']:+.6f}); "
          f"row-matched CD4 macro {100*S['arms'][f'{arm}CD4ROWS@0.35']['rel_macro']:+.3f}% (agg "
          f"{100*S['arms'][f'{arm}CD4ROWS@0.35']['rel_agg']:+.3f}%, dpds {S['arms'][f'{arm}CD4ROWS@0.35']['d_pds']:+.6f})")
        rs = a["rowsplit"]
        P(f"  rowsplit: carrying rows rel cos {100*rs['rel_in']:+.3f}%, NOT-carrying {100*rs['rel_out']:+.5f}% "
          f"(largest |excursion| {rs['max_abs_rel_out_pp']:.5f} pp vs 218's 0.018 pp bar)")
        P(f"  leave-top-3-out: pds retained {100*conc['pds_retained']:.1f}% (top-3 share {100*conc['top_share_pds']:.1f}%), "
          f"cos retained {100*conc['cos_retained']:.1f}%")
        P(f"  delivered amplitude vs h1union: r = {r_amp:.5f} ({100*(r_amp-1):+.3f}%) from {amp_src}")
        P(f"  216 identity: nmae wald89 {six['nmae_wald89']:.5f} (incumbent {SIX['h1union g=0.28 <- INCUMBENT']['nmae_wald89']:.5f}), "
          f"d {six['d_nmae_vs_incumbent']:+.5f} raw -> overall {six['overall_nmae_vs_incumbent']:+.6f}; "
          f"e = dM/dU {six['efficiency_dM_over_dU']:+.3f} (break-even 0.5); d top-89 precision {six['d_precision_vs_incumbent']:+.5f}")

        rows = {}
        P(f"\n  (i) 243's class rule -- direction legs + DE-four propagation:")
        P(f"    {'class':52s}{'k_cos':>6s}{'k_pds':>6s}{'d cos':>10s}{'raw mse':>9s}{'o_mse':>9s}{'o_pds':>9s}{'o_DE4':>9s}"
          f"{'TOTAL (i)':>11s}{'x gate':>8s}{'nmae id':>9s}{'TOTAL (ii)':>11s}")
        for label, (kc, kp, note) in CLASSES.items():
            dcos = kc * rel_agg * cos_inc
            cos_arm = cos_inc + dcos
            if per_arm is not None:
                raw_arm = raw_mse(per_arm, cos_arm)
            else:
                raw_arm = raw_mse(per_inc, cos_arm, r=r_amp)
            o_mse = -(raw_arm - raw_inc) / SPAN_MSE / 6.0
            o_pds = kp * d_pds / SPAN_PDS / 6.0
            o_de4 = K_DE4["h1union point"] * dcos
            total_i = o_mse + o_pds + o_de4
            o_nmae_id = six["overall_nmae_vs_incumbent"]
            total_ii = o_mse + o_pds + o_de4 * (1 - NMAE_SHARE_OF_DE4) + o_nmae_id
            rows[label] = dict(k_cos=kc, k_pds=kp, d_cos=dcos, cos_arm=cos_arm, raw_mse=raw_arm,
                               o_mse=o_mse, o_pds=o_pds, o_de4=o_de4, total_i=total_i,
                               o_nmae_identity=o_nmae_id, total_ii=total_ii, x_gate=total_i / GATE)
            P(f"    {label:52s}{kc:>6.2f}{kp:>6.2f}{dcos:>+10.6f}{raw_arm:>9.5f}{o_mse:>+9.5f}{o_pds:>+9.5f}{o_de4:>+9.5f}"
              f"{total_i:>+11.5f}{total_i/GATE:>8.2f}{o_nmae_id:>+9.5f}{total_ii:>+11.5f}")
        # the amplitude leg alone (k_cos = 0), for the record
        raw_amp_only = raw_mse(per_arm, cos_inc) if per_arm is not None else raw_mse(per_inc, cos_inc, r=r_amp)
        P(f"    amplitude/budget leg alone (k_cos 0): raw mse {raw_amp_only:.5f} (d {raw_amp_only-raw_inc:+.5f}), "
          f"overall {-(raw_amp_only-raw_inc)/SPAN_MSE/6:+.6f}")
        # DE-four sensitivity at the point class
        pt = rows["new-source class, point (243 §4.4)"]
        P(f"    DE-four sensitivity (point class): low rate {K_DE4['low (a14->a15)']*pt['d_cos']:+.5f}, "
          f"high rate {K_DE4['high (cd4f034)']*pt['d_cos']:+.5f} (carried {pt['o_de4']:+.5f})")

        # the row-matched CD4 control and the two shuffles under the SAME rule (point class), so the
        # falsifiers are priced in the units the gate reads, not only in PLOLO units
        ctl = {}
        for vtag in ("CD4ROWS", "SHUF20726", "SHUF771103"):
            c = S["arms"][f"{arm}{vtag}@0.35"]
            kc, kp, _n = CLASSES["new-source class, point (243 §4.4)"]
            dcos_c = kc * c["rel_agg"] * cos_inc
            raw_c = raw_mse(per_inc, cos_inc + dcos_c, r=r_amp)
            o_mse_c = -(raw_c - raw_inc) / SPAN_MSE / 6.0
            o_pds_c = kp * c["d_pds"] / SPAN_PDS / 6.0
            o_de4_c = K_DE4["h1union point"] * dcos_c
            ctl[vtag] = dict(rel_agg=c["rel_agg"], rel_macro=c["rel_macro"], d_pds=c["d_pds"], d_cos=dcos_c,
                             o_mse=o_mse_c, o_pds=o_pds_c, o_de4=o_de4_c, total_i=o_mse_c + o_pds_c + o_de4_c)
            P(f"    control {vtag:11s} under the point class: rel agg {100*c['rel_agg']:+.3f}% macro {100*c['rel_macro']:+.3f}% "
              f"dpds {c['d_pds']:+.6f} -> o_mse {o_mse_c:+.5f} o_pds {o_pds_c:+.5f} o_DE4 {o_de4_c:+.5f} TOTAL {ctl[vtag]['total_i']:+.5f}")
        point = rows["new-source class, point (243 §4.4)"]
        P(f"    the real source beats the row-matched CD4 control on the composite by "
          f"{point['total_i'] - ctl['CD4ROWS']['total_i']:+.5f} (on cos_agg alone the control is {'LARGER' if ctl['CD4ROWS']['rel_agg'] > rel_agg else 'smaller'}: "
          f"{100*ctl['CD4ROWS']['rel_agg']:+.3f}% vs {100*rel_agg:+.3f}%; on macro cos and pds the source wins)")
        gate_i = point["total_i"] >= GATE_X * GATE
        gate_pds = a["pds_pos"] == 3 and d_pds >= 0
        gate_conc = conc["pds_retained"] >= 0.5
        verdict = gate_i and gate_pds and gate_conc
        P(f"\n  GATE (mandate C) for {arm}:")
        P(f"    (i) point class {point['total_i']:+.5f} >= 1.5 x 0.0019 = {GATE_X*GATE:.5f}: {'PASS' if gate_i else 'FAIL'} "
          f"({point['total_i']/GATE:.2f}x the shape-only gate; best class {max(r['total_i'] for r in rows.values()):+.5f})")
        P(f"    d pds >= 0 on 3/3 cancer folds: {'PASS' if gate_pds else 'FAIL'} ({a['pds_pos']}/3, mean {d_pds:+.6f})")
        P(f"    leave-top-3-out pds retention >= 50%: {'PASS' if gate_conc else 'FAIL'} ({100*conc['pds_retained']:.1f}%)")
        P(f"    => {'CLEARS -- build' if verdict else 'DOES NOT CLEAR -- do not build'}")
        results[arm] = dict(plolo=dict(rel_agg=rel_agg, rel_macro=rel_macro, d_pds=d_pds, pds_pos=a["pds_pos"],
                                       d_qw=a["d_qw"], per_fold=a["per_fold"], kolf2=k2),
                            amplitude=dict(r=r_amp, source=amp_src),
                            nmae_identity=dict(d_raw=six["d_nmae_vs_incumbent"], overall=six["overall_nmae_vs_incumbent"],
                                               e=six["efficiency_dM_over_dU"], d_prec=six["d_precision_vs_incumbent"]),
                            conc3=conc, rows=rows, controls_point_class=ctl, gate=dict(i=gate_i, pds=gate_pds, conc=gate_conc, clears=verdict,
                                                              point_total_i=point["total_i"], threshold=GATE_X * GATE))
    (WORK / "price248.json").write_text(json.dumps(dict(results=results, cos_inc=cos_inc, raw_inc=raw_inc,
                                                         classes=CLASSES, k_de4=K_DE4,
                                                         nmae_share_of_de4=NMAE_SHARE_OF_DE4), indent=1, default=float))
    (WORK / "price248.log").write_text("\n".join(lines) + "\n")
    P(f"\nwrote {WORK/'price248.json'}")


if __name__ == "__main__":
    main()
