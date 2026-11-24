"""254 -- the six-member price of the composed arm and of its INCREMENT over cd4svd5, every member
explicit, never at zero (246 section 4's route = 218 section 7.4's rate card; 243's class rule as
the sensitivity).

Route (246 `price246.py`, verbatim constants from `price218`):
  * the five cosine-driven members (`mse`, `reach`, `nmae`, `fid`, `jac`) from `cd4f034`'s realised
    board rate card per +0.01 ABSOLUTE cosine, on the MACRO PLOLO cosine delta at k_cos = 1;
  * `pds` from PLOLO directly (raw d / 0.452 / 6 x 1.011) at THREE transfer factors 0.7x / 1.0x / 2.0x.
KNOWN-ANSWER ROWS: 246's own `E_cd4_svd5+partial@0.34`, `B2_cd4_svd5_resolved` and `C2_partial@0.34`
prices (price246.json) reproduce from this experiment's fold numbers to 1e-9.

    uv run python price254.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
from price218 import BOARD_COS, GATE_AMPL, GATE_SHAPE, RATE, SIGMA_PDS, SPAN_PDS  # noqa: E402

E246 = ROOT / "experiments/246-direction-hypotheses"
NOISE_ORDERING = 0.004502     # AGENTS.md rule 13: sd of board d overall between digest-identical arms


def six(rel_macro: float, d_pds: float) -> dict:
    dcos = BOARD_COS * rel_macro
    legs = {m: RATE[m] * dcos / 0.01 for m in RATE if m != "pds"}
    pds1 = d_pds / SPAN_PDS / 6.0 * 1.011
    tot = {f"{f:g}": sum(legs.values()) + f * pds1 for f in (0.7, 1.0, 2.0)}
    return dict(rel_macro=rel_macro, d_pds=d_pds, pds_sigma=abs(d_pds) / SIGMA_PDS, d_cos_abs=dcos,
                legs=legs, pds={f"{f:g}": f * pds1 for f in (0.7, 1.0, 2.0)}, totals=tot)


def main() -> None:
    S = json.loads((WORK / "repro254-summary.json").read_text())
    P246 = json.loads((E246 / "work/price246.json").read_text())
    six246 = json.loads((E246 / "work/six246.json").read_text())["arms"]
    lines = []

    def P(t=""):
        print(t)
        lines.append(t)

    arms = {
        "E_cd4_svd5+partial@0.34 (246's construction, FIXED fraction; KA)": S["composed_fixed"],
        "B2_cd4_svd5_resolved = cd4svd5 (249's arm; KA)": S["svd5"],
        "C2_partial@0.34 alone (KA)": S["partial"],
        "COMPOSED resolved (the package) vs h1union": S["composed_resolved"],
        "INCREMENT: composed vs cd4svd5 (the slot's one variable)": S["increment_over_svd5"],
    }
    ka_map = {
        "E_cd4_svd5+partial@0.34 (246's construction, FIXED fraction; KA)": "E_cd4_svd5+partial@0.34",
        "B2_cd4_svd5_resolved = cd4svd5 (249's arm; KA)": "B2_cd4_svd5_resolved",
        "C2_partial@0.34 alone (KA)": "C2_partial@0.34",
    }
    out = {}
    P(f"{'arm':66s}{'macro%':>8s}{'d pds':>10s}{'sig':>5s}{'mse':>9s}{'reach':>9s}{'nmae':>9s}{'fid':>9s}"
      f"{'jac':>9s}{'pds x0.7':>9s}{'pds x1':>8s}{'pds x2':>8s}{'TOT x0.7':>9s}{'TOT x1':>8s}{'TOT x2':>8s}")
    ka_ok = True
    for label, s in arms.items():
        r = six(s["rel_macro"], s["d_pds"])
        out[label] = r
        P(f"{label:66s}{100*r['rel_macro']:+8.3f}{r['d_pds']:+10.6f}{r['pds_sigma']:5.1f}"
          f"{r['legs']['mse']:+9.5f}{r['legs']['reach']:+9.5f}{r['legs']['nmae']:+9.5f}{r['legs']['fid']:+9.5f}"
          f"{r['legs']['jac']:+9.5f}{r['pds']['0.7']:+9.5f}{r['pds']['1']:+8.5f}{r['pds']['2']:+8.5f}"
          f"{r['totals']['0.7']:+9.5f}{r['totals']['1']:+8.5f}{r['totals']['2']:+8.5f}")
        if label in ka_map:
            ref = P246[ka_map[label]]
            d = max(abs(r["totals"]["0.7"] - ref["totals"]["0.7"]), abs(r["totals"]["1"] - ref["totals"]["1.0"]),
                    abs(r["totals"]["2"] - ref["totals"]["2.0"]))
            ok = d < 1e-9
            ka_ok &= ok
            r["ka_vs_price246"] = dict(ref_totals=ref["totals"], max_abs_diff=d, ok=ok)
            P(f"{'':66s}  KA vs price246 {ka_map[label]}: totals {ref['totals']['0.7']:+.5f} / {ref['totals']['1.0']:+.5f} / "
              f"{ref['totals']['2.0']:+.5f}  max|d| {d:.1e} {'PASS' if ok else 'FAIL'}")
    P(f"\ngates: shape-only {GATE_SHAPE:+.4f}   amplitude-changing {GATE_AMPL:+.4f}   PLOLO pds LOO sigma {SIGMA_PDS} raw   "
      f"ordering noise sd between digest-identical arms {NOISE_ORDERING} overall")
    inc = out["INCREMENT: composed vs cd4svd5 (the slot's one variable)"]
    comp = out["COMPOSED resolved (the package) vs h1union"]
    P(f"\nthe INCREMENT against the shape-only gate: {inc['totals']['0.7']/GATE_SHAPE:.2f}x / {inc['totals']['1']/GATE_SHAPE:.2f}x / "
      f"{inc['totals']['2']/GATE_SHAPE:.2f}x at pds transfer 0.7 / 1.0 / 2.0  -- BELOW the gate at every factor "
      f"(246 section 3.7: NOT build-worthy alone); leave-top-3-out of the increment "
      f"{100*S['concentration']['increment']['pds_retained']:+.1f}% (bar 50%)")
    P(f"the COMPOSED arm against h1union: {comp['totals']['0.7']/GATE_SHAPE:.2f}x / {comp['totals']['1']/GATE_SHAPE:.2f}x / "
      f"{comp['totals']['2']/GATE_SHAPE:.2f}x the shape-only gate; leave-top-3-out "
      f"{100*S['concentration']['composed_resolved']['pds_retained']:.1f}%")
    ex_e = six246.get("E_cd4_svd5+partial@0.34", {})
    ex_c = six246.get("C_partial@0.34", {})
    P(f"\nexact `nmae` identity (six246, a magnitude bound, not a sign): E {ex_e.get('overall_nmae_vs_inc', float('nan')):+.5f} "
      f"(e {ex_e.get('e_vs_inc', float('nan')):.3f}), C2 alone {ex_c.get('overall_nmae_vs_inc', float('nan')):+.5f} "
      f"(e {ex_c.get('e_vs_inc', float('nan')):.3f}); delivered ratio vs a15: E {ex_e.get('delivered_norm_ratio')}, C2 {ex_c.get('delivered_norm_ratio')}")
    # board-raw translation of the pds legs (span 0.452, overall = raw / 0.452 / 6 x 1.011 -> raw = d_pds x factor)
    P(f"\nboard raw pds implied for the INCREMENT at 0.7 / 1.0 / 2.0x: {0.7*inc['d_pds']:+.5f} / {inc['d_pds']:+.5f} / {2*inc['d_pds']:+.5f} "
      f"on top of cd4svd5's; for the COMPOSED arm vs h1union's 0.83572: {0.7*comp['d_pds']:+.5f} / {comp['d_pds']:+.5f} / {2*comp['d_pds']:+.5f}")
    (WORK / "price254.json").write_text(json.dumps(dict(arms=out, ka_all=bool(ka_ok), gates=dict(shape=GATE_SHAPE, amplitude=GATE_AMPL),
                                                        noise_ordering=NOISE_ORDERING), indent=1))
    (WORK / "price254.log").write_text("\n".join(lines) + "\n")
    print(f"\nKA {'PASS' if ka_ok else 'FAIL'}; wrote {WORK/'price254.json'}")
    raise SystemExit(0 if ka_ok else 1)


if __name__ == "__main__":
    main()
