"""215 Phase 0 — the five known-answer arms.  ABORT if the decomposition misclassifies any.

Every number below is a BOARD number with its provenance named.  Nothing is modelled.
"""
import json, os, sys
import numpy as np
import chan215 as C

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
os.makedirs(WORK, exist_ok=True)
M = C.MEMBERS


def from_scaled_delta(anchor_raw, dsca):
    """raw of a sibling arm = anchor raw + delta(scaled)/slope.  Cancels the intercept error."""
    return {m: anchor_raw[m] + dsca[m] / C.SLOPE[m] for m in M}


# ---------------------------------------------------------------------------- board rows
# a14, a15: docs/BOARD-LEDGER.md "Raw members, where recorded"  (all six recorded)
A14 = dict(pds=0.832702, mse=0.882775, nmae=0.945764, fid=0.502881, reach=0.180615, jac=0.030081)
A15 = dict(pds=0.835987, mse=0.873805, nmae=0.944251, fid=0.503536, reach=0.181173, jac=0.030080)

# rc_hybrid, rc_callset: the LIVE board (work/lb-live.json), full float precision
LB = json.load(open(os.path.join(WORK, "lb-live.json")))["entries"]
RAWK = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
        "fid": "de_wilcoxon_direction_fidelity_yield_raw",
        "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard"}
def board_row(model):
    e = [x for x in LB if x["model_name"] == model and x["team_name"] == "Ahmet Dedeler"]
    assert len(e) == 1, f"{model}: {len(e)} rows"
    return {m: float(e[0][RAWK[m]]) for m in M}
RCC = board_row("Peer Review Rejected This")            # rc_callset, entry 31oB9yWy8CSNBDEkY4iW
# rc_hybrid is no longer the displayed row, so it is not on the live board; take it from the panel
PANEL = json.load(open(os.path.join(WORK, "panel215.json")))
_rh = [r for r in PANEL if r["team"] == "Ahmet Dedeler"
       and r["model"] == "The Mitochondria Declined To Comment"][0]
RCH = {m: _rh["raw_" + m] for m in M}

# a11: BOARD-LEDGER raw row (pds/mse/reach) + 197 §1 raw mse; a10 built from the scaled delta
A11 = dict(pds=0.833307, mse=0.909989,
           nmae=1.00034 + 0.096982 / C.SLOPE["nmae"],
           fid=0.51241 + (-0.031535) / C.SLOPE["fid"],
           reach=0.178243,
           jac=0.03037 + (-0.002907) / C.SLOPE["jac"])
# a10 -> a11 scaled deltas, docs/BOARD-LEDGER.md
D_A10_A11 = dict(pds=0.739889 - 0.736561, mse=None, nmae=0.096982 - 0.094338,
                 fid=-0.031535 - (-0.027897), reach=0.112220 - 0.123595,
                 jac=-0.002907 - (-0.001493))
A10 = {m: A11[m] - D_A10_A11[m] / C.SLOPE[m] for m in M if m != "mse"}
A10["mse"] = 1.352415        # 197 §1, board raw (a11's 0.909989 is in the ledger too)

# a19: docs/LOG.md 2026-08-31 BOARD RESULT table, scaled deltas against a15
D_A15_A19 = dict(pds=-0.00118, mse=-0.09591, nmae=+0.00845,
                 fid=+0.00400, reach=+0.00297, jac=+0.00057)
A19 = from_scaled_delta(A15, D_A15_A19)

ARMS = [
    ("a10 -> a11        alpha 1.0796 -> 0.838", A10, A11, "AMPLITUDE",
     "197 §1: cosine 0.07554 -> 0.07550, i.e. FLAT"),
    ("a14 -> a15        CD4 residual",          A14, A15, "DIRECTION",
     "197 §1.2: cosine 0.07418 -> 0.07925, +6.8%"),
    ("a15 -> a19        whitener ridge",        A15, A19, "GEOMETRY",
     "LOG 08-31: mse -0.09591 scaled, pds -0.00118, nmae IMPROVED"),
    ("a15 -> rc_hybrid  rank-controlled emit",  A15, RCH, "ORDERING",
     "pds and mse reproduce a15 to 16 s.f."),
    ("rc_hybrid -> rc_callset  call set",       RCH, RCC, "ORDERING",
     "pds and mse bit-identical to rc_hybrid"),
]


def main():
    rows, ok = [], True
    print("=" * 108)
    print("PHASE 0 — the five known-answer arms")
    print("=" * 108)
    hdr = f"{'arm':<40}{'dpds':>10}{'dmse':>10}{'dnmae':>10}{'dreach':>10}{'rho':>9}{'dcos_pds':>10}{'dcos_mse':>10}"
    print(hdr)
    for name, b, a, expect, note in ARMS:
        r = C.classify(b, a)
        rows.append({"arm": name, "expect": expect, "note": note, **{
            k: r[k] for k in ("label", "rho", "dcos_pds", "dcos_mse", "geom_share",
                              "b_direction", "b_ordering", "b_amplitude", "b_geometry",
                              "ovr_obs", "both_informative", "dcos_pds_sd", "sigma_pds_eff")}})
        print(f"{name:<40}{r['d_raw']['pds']:>+10.5f}{r['d_raw']['mse']:>+10.5f}"
              f"{r['d_raw']['nmae']:>+10.5f}{r['d_raw']['reach']:>+10.5f}"
              f"{r['rho']:>9.3f}{r['dcos_pds']:>+10.5f}{r['dcos_mse']:>+10.5f}")

    print("\nCONFUSION MATRIX  (rows = truth, cols = what the instrument said)")
    labels = ["DIRECTION", "DIRECTION+AMPLITUDE", "ORDERING", "AMPLITUDE", "GEOMETRY", "NULL"]
    print(f"{'truth':<12}" + "".join(f"{l[:11]:>13}" for l in labels))
    truth_set = ["DIRECTION", "ORDERING", "AMPLITUDE", "GEOMETRY"]
    for t in truth_set:
        counts = {l: 0 for l in labels}
        for rr in rows:
            if rr["expect"] == t:
                counts[rr["label"]] += 1
        print(f"{t:<12}" + "".join(f"{counts[l]:>13d}" for l in labels))
    for rr in rows:
        good = rr["label"].startswith(rr["expect"]) or rr["expect"] in rr["label"]
        ok &= good
        rr["pass"] = bool(good)

    print("\nPER-ARM VERDICT")
    print(f"{'arm':<40}{'truth':<11}{'called':<21}{'verdict':<7}  known answer")
    for rr in rows:
        print(f"{rr['arm']:<40}{rr['expect']:<11}{rr['label']:<21}"
              f"{'PASS' if rr['pass'] else 'FAIL':<7}  {rr['note']}")

    print("\nCHANNEL BUDGET (overall units; the four buckets sum EXACTLY to the observed delta)")
    print(f"{'arm':<40}{'direction':>11}{'ordering':>11}{'amplitude':>11}{'geometry':>11}{'sum':>11}{'observed':>11}")
    for rr in rows:
        s = rr["b_direction"] + rr["b_ordering"] + rr["b_amplitude"] + rr["b_geometry"]
        print(f"{rr['arm']:<40}{rr['b_direction']:>+11.5f}{rr['b_ordering']:>+11.5f}"
              f"{rr['b_amplitude']:>+11.5f}{rr['b_geometry']:>+11.5f}{s:>+11.5f}{rr['ovr_obs']:>+11.5f}")
        assert abs(s - rr["ovr_obs"]) < 1e-12, "budget does not partition the observed delta"

    print("\nMAGNITUDE CHECK — the direction channel against 197's independently inverted cosine")
    mag = [
        ("a10 -> a11", rows[0]["dcos_pds"], -0.00004, "197 §1: 0.07554 -> 0.07550"),
        ("a14 -> a15", rows[1]["dcos_pds"], +0.00507, "197 §1.2: 0.07418 -> 0.07925"),
        ("a15 -> rc_hybrid", rows[3]["dcos_pds"], +0.000001, "197 §5 falsifier 4: 2.5e-7"),
        ("rc_hybrid -> rc_callset", rows[4]["dcos_pds"], 0.0, "202: pds bit-identical"),
    ]
    for (nm, got, want, src), rr in zip(mag, [rows[0], rows[1], rows[3], rows[4]]):
        tol = rr.get("dcos_pds_sd", 0.00232)
        z = abs(got - want) / tol
        good = z <= 1.5
        print(f"  {nm:<26} instrument {got:>+9.5f}   197 {want:>+9.5f}   "
              f"|err| {abs(got-want):.5f} = {z:4.2f} sigma (floor {tol:.5f})   "
              f"{'PASS' if good else 'FAIL'}   {src}")
        ok &= good

    print()
    if not ok:
        print("PHASE 0 FAILED — the instrument is not calibrated.  Nothing downstream is readable.")
        sys.exit(1)
    print("PHASE 0: 5 of 5 arms correctly classified, 4 of 4 magnitudes inside the pds floor.")
    json.dump(rows, open(os.path.join(WORK, "phase0215.json"), "w"), indent=1, default=float)
    # export the arm table for RESULTS.md
    json.dump({"A10": A10, "A11": A11, "A14": A14, "A15": A15, "A19": A19,
               "RCH": RCH, "RCC": RCC, "E_DIR": C.E_DIR, "RATIO_DIR": C.RATIO_DIR,
               "RATIO_DIR_MSE_IDENTITY": C.RATIO_DIR_MSE_IDENTITY, "SIGMA": C.SIGMA},
              open(os.path.join(WORK, "calib215.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
