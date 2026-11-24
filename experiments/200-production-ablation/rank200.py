"""200 — assemble the six-member leave-one-out table in OVERALL board units, and rank.

Sign convention throughout: **positive = the component EARNS its place**, i.e. the number is
what we would LOSE by removing it.  For a member measured as an ablation delta (arm minus
base) the sign is flipped here once, in `LOO`, and nowhere else.

The scoring map is the board-fitted one (`docs/BOARD-LEDGER.md`, 612 live entries):
`dOverall/dRaw` = pds +0.3690, mse -0.1752, nmae -0.2755, fid +0.5607, reach +0.1882,
jac +0.4544.  The per-member numbers below arrive already in overall units from their own
instrument, each of which applies that map (`report_pds.py`, `mse200.py`, `six200.py`).

Gate: rule 13's composed slot gate, +0.0019 for a shape-only arm and +0.0089 for an
amplitude-changing one.  A component whose removal scores BELOW -(gate) is a candidate for
removal; one whose removal scores above +(gate) is load-bearing; anything between is inside
the unpriced-channel noise and is called noise, not value.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")

GATE_SHAPE, GATE_AMP = 0.0019, 0.0089

# component -> (pds key in report_pds.json loo, mse source, DE key in six200-table.json)
# `mse` sources: ("emission", <label>) reads mse200.json["emission"], which is the VALIDATED
# route; ("geometry", <arm>) reads mse200.json["geometry"], which is the rX = 1 partial price.
COMPONENTS = [
    ("union whitening operator (gamma 0.65)", "noWhiten", ("geometry", "noWhiten"),
     "whitening (gamma 0.65)", "shape"),
    # DE key is None ON PURPOSE.  h1-test has CD4 coverage on 4 of 85 targets, and production's
    # `solve_cd4_fraction` sizes the residual by the MEDIAN row norm of the added block -- with
    # 4/85 rows non-zero that median is exactly 0, so the solved fraction is 0.0 and the arm is
    # BYTE-IDENTICAL to `reb_ship` (verified: sha256 aa28e8495ab66101... on all three seeds).
    # A structural zero is not a measurement, so CD4 is priced on `pds` + `mse` only.
    ("post-whitening CD4 residual (f = 0.10)", "noCD4", ("geometry", "noCD4"),
     None, "amp"),
    ("PerturbMatch K562 replacement (0.5)", "noPM", ("geometry", "noPM"),
     "PerturbMatch blend 0.5", "shape"),
    ("equal donor weights (1,1,1)", "weights_4_1_1", ("geometry", "weights_4_1_1"),
     "equal donor weights", "shape"),
    ("CPM gate (--min-cpm 10)", "noGate", ("geometry", "noGate"),
     "CPM gate 10 (shipped) LOO", "shape"),
    ("pruned basis (union_pruned)", "noPrune", ("geometry", "noPrune"), None, "shape"),
    ("VCC_ALPHA 0.982 (vs v11's 1.0796)", "alpha_v11", ("geometry", "alpha_v11"),
     None, "amp"),
    ("a23 mass fix (minnorm)", "noMassFix", ("geometry", "noMassFix"), "a23 mass fix", "shape"),
    ("tiered spread policy (kt 8000/1.0)", None, ("emission", "-tier (flat 0.25) phi4"),
     "tiered spread policy", "shape"),
    ("--overdispersion 4.0", None, ("emission", "-phi  tier8000 phi0"),
     "--overdispersion 4.0", "shape"),
]


# The same arithmetic for knobs that were SET rather than added -- the shipped value against
# its nearest measured neighbours.  Sign convention identical: positive = the SHIPPED setting
# is worth that much more than the alternative.
KNOBS = [
    ("WHITEN_GAMMA 0.65 vs 0.80", "gamma_080", ("geometry", "gamma_080"),
     "gamma 0.65 vs 0.80", "shape"),
    ("WHITEN_GAMMA 0.65 vs 0.50", "gamma_050", ("geometry", "gamma_050"),
     "gamma 0.65 vs 0.50", "shape"),
    ("donor weights 1:1:1 vs 2:1:1", "weights_2_1_1", ("geometry", "weights_2_1_1"),
     "equal vs 2:1 K562-heavy", "shape"),
    ("CPM gate 10 vs 5 (088's move)", "gate5", ("geometry", "gate5"),
     "CPM gate 10 vs 5 (088)", "shape"),
    ("CPM gate 10 vs 20", "gate20", ("geometry", "gate20"),
     "CPM gate 10 vs 20", "shape"),
    ("tier 8000 vs FLAT 1.00", None, ("emission", "flat 1.00 phi4"),
     "tier8000 vs flat 1.00", "shape"),
    ("tier 8000 vs tier 4000", None, ("emission", "tier 4000 phi4"),
     "tier8000 vs tier4000", "shape"),
    ("tier 8000 vs tier 14000", None, ("emission", "tier 14000 phi4"), None, "shape"),
    ("phi 4 vs phi 8", None, ("emission", "tier8000 phi 8"), None, "shape"),
]


def main() -> None:
    pds = json.loads(open(os.path.join(WORK, "report_pds.json")).read())
    mse = json.loads(open(os.path.join(WORK, "mse200.json")).read())
    six = {}
    p = os.path.join(WORK, "six200-tablev2.json")
    if os.path.exists(p):
        six = json.loads(open(p).read())

    rows = []
    for label, pkey, (msrc, mkey), dekey, cls in COMPONENTS + KNOBS:
        r = dict(component=label, cls=cls)
        # ---- pds
        if pkey is None:
            r["pds"] = 0.0
            r["pds_note"] = "0 by construction (emission-only; board-verified twice)"
        elif pkey in pds["loo"]:
            r["pds"] = -pds["loo"][pkey]["overall_k1011"]
            r["pds_sd"] = pds["loo"][pkey]["sd"] / 0.45208 * 1.011 / 6
            r["pds_note"] = f"PLOLO 4-fold, signs {pds['loo'][pkey]['signs']}"
        else:
            r["pds"] = None
        # ---- mse
        if msrc == "emission" and mkey in mse.get("emission", {}):
            r["mse"] = -mse["emission"][mkey]["d_overall"]
            r["mse_note"] = ("096 claim route, VALIDATED class "
                             f"(claim {mse['emission'][mkey]['claim']:.5f} vs "
                             f"s^2B 0.23268, branch {mse['emission'][mkey]['branch']})")
        elif msrc == "geometry" and mkey in mse.get("geometry", {}):
            g = mse["geometry"][mkey]
            r["mse"] = -g["d_overall"]
            r["mse_note"] = (f"partial: rX = 1 null, rA {g['rA']:.5f} rB {g['rB']:.5f} "
                             f"rClaim {g['rClaim']:.5f}")
        else:
            r["mse"] = None
        # ---- the four DE members
        if dekey and dekey in six:
            d = six[dekey]
            for m in ("reach", "fid", "jac", "nmae"):
                r[m] = d[m]
            r["jac_np"] = d.get("jac_np")
            r["de_note"] = (f"h1-test, {d['paired_reach']['n']} matched seeds, "
                            f"d reach raw {d['d_reach_raw']:+.5f} "
                            f"(sem {d['paired_reach']['sem']:.5f})")
        else:
            for m in ("reach", "fid", "jac", "nmae"):
                r[m] = None
            r["jac_np"] = None
        known = [r[m] for m in ("pds", "mse", "reach", "fid", "jac", "nmae")
                 if r[m] is not None]
        r["sum_known"] = sum(known)
        r["n_known"] = len(known)
        # the same sum with the call-set-scaled `jac` convention (six200.overall docstring)
        r["sum_jacnp"] = (r["sum_known"] - (r["jac"] or 0.0) + (r["jac_np"] or 0.0))
        r["knob"] = label in {k[0] for k in KNOBS}
        rows.append(r)

    rows.sort(key=lambda x: (x["knob"], -x["sum_known"]))
    W = 40
    print("=" * 132)
    print("SIX-MEMBER LEAVE-ONE-OUT   (positive = what removing the component would COST us, "
          "in overall board units)")
    print("=" * 132)
    print(f"{'component':<{W}}{'pds':>10}{'mse':>10}{'reach':>10}{'fid':>10}{'jac':>10}"
          f"{'nmae':>10}{'SUM':>11}{'SUM jac(np)':>13}{'n':>4}  class")
    fmt = lambda v: "     --  " if v is None else f"{v:>+10.6f}"  # noqa: E731
    seen_knob = False
    for r in rows:
        if r["knob"] and not seen_knob:
            seen_knob = True
            print("-" * 132)
            print("KNOB NEIGHBOURHOOD -- the shipped setting against its nearest "
                  "measured alternative (not a leave-one-out)")
            print("-" * 132)
        print(f"{r['component']:<{W}}" + "".join(fmt(r[m]) for m in
              ("pds", "mse", "reach", "fid", "jac", "nmae"))
              + f"{r['sum_known']:>+11.6f}{r['sum_jacnp']:>+13.6f}{r['n_known']:>4}  "
              f"{r['cls']}")

    print()
    print("=" * 132)
    print("VERDICT   (gate: shape-only +/-0.0019, amplitude-changing +/-0.0089)")
    print("=" * 132)
    for r in rows:
        if r["knob"]:
            continue
        gate = GATE_SHAPE if r["cls"] == "shape" else GATE_AMP
        s = r["sum_known"]
        s2 = r["sum_jacnp"]
        v = ("LOAD-BEARING" if s > gate else
             "COSTING US -- REMOVE" if s < -gate else "NOISE (inside the gate)")
        v2 = ("LOAD-BEARING" if s2 > gate else
              "COSTING US -- REMOVE" if s2 < -gate else "NOISE (inside the gate)")
        flag = "" if v == v2 else f"   <-- FLIPS under jac(np): {v2}"
        print(f"  {r['component']:<{W}}{s:>+11.6f}   vs +/-{gate:.4f}   {v}{flag}")

    print()
    print("notes per member, so nothing above is read without its instrument:")
    for r in rows:
        print(f"  {r['component']}")
        for k in ("pds_note", "mse_note", "de_note"):
            if r.get(k):
                print(f"      {k[:-5]:<6} {r[k]}")

    json.dump(rows, open(os.path.join(WORK, "rank200.json"), "w"), indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'rank200.json')}")


if __name__ == "__main__":
    sys.exit(main())
