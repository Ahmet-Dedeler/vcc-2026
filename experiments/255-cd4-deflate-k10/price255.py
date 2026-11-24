"""255 -- the six-member price of `cd4svd10` vs `h1union` and of its INCREMENT over `cd4svd5` (the
slot's one variable), every member explicit, never at zero (246 section 4's route = 218 section
7.4's rate card; 254's `price254.py` form).

Route (246 `price246.py`, verbatim constants from `price218`):
  * the five cosine-driven members (`mse`, `reach`, `nmae`, `fid`, `jac`) from `cd4f034`'s realised
    board rate card per +0.01 ABSOLUTE cosine, on the MACRO PLOLO cosine delta at k_cos = 1;
  * `pds` from PLOLO directly (raw d / 0.452 / 6 x 1.011) at THREE transfer factors 0.7x / 1.0x / 2.0x.
KNOWN-ANSWER ROWS: 246's own `B2_cd4_svd10_resolved` and `B2_cd4_svd5_resolved` prices
(price246.json) reproduce from this experiment's fold numbers to 1e-9.

    uv run python price255.py
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
    S = json.loads((WORK / "repro255-summary.json").read_text())
    P246 = json.loads((E246 / "work/price246.json").read_text())
    six246 = json.loads((E246 / "work/six246.json").read_text())["arms"]
    lines = []

    def P(t=""):
        print(t)
        lines.append(t)

    arms = {
        "cd4svd10 = B2_cd4_svd10_resolved vs h1union (KA)": S["arm"],
        "cd4svd5  = B2_cd4_svd5_resolved vs h1union (KA)": S["arm5"],
        "INCREMENT: cd4svd10 vs cd4svd5 (the slot's one variable)": S["increment"],
    }
    ka_map = {
        "cd4svd10 = B2_cd4_svd10_resolved vs h1union (KA)": "B2_cd4_svd10_resolved",
        "cd4svd5  = B2_cd4_svd5_resolved vs h1union (KA)": "B2_cd4_svd5_resolved",
    }
    out = {}
    P(f"{'arm':60s}{'macro%':>8s}{'d pds':>10s}{'sig':>5s}{'mse':>9s}{'reach':>9s}{'nmae':>9s}{'fid':>9s}"
      f"{'jac':>9s}{'pds x0.7':>9s}{'pds x1':>8s}{'pds x2':>8s}{'TOT x0.7':>9s}{'TOT x1':>8s}{'TOT x2':>8s}")
    ka_ok = True
    for label, s in arms.items():
        r = six(s["rel_macro"], s["d_pds"])
        out[label] = r
        P(f"{label:60s}{100*r['rel_macro']:+8.3f}{r['d_pds']:+10.6f}{r['pds_sigma']:5.1f}"
          f"{r['legs']['mse']:+9.5f}{r['legs']['reach']:+9.5f}{r['legs']['nmae']:+9.5f}{r['legs']['fid']:+9.5f}"
          f"{r['legs']['jac']:+9.5f}{r['pds']['0.7']:+9.5f}{r['pds']['1']:+8.5f}{r['pds']['2']:+8.5f}"
          f"{r['totals']['0.7']:+9.5f}{r['totals']['1']:+8.5f}{r['totals']['2']:+8.5f}")
        if label in ka_map:
            ref = P246[ka_map[label]]
            d = max(abs(r["totals"]["0.7"] - ref["totals"]["0.7"]), abs(r["totals"]["1"] - ref["totals"]["1.0"]),
                    abs(r["totals"]["2"] - ref["totals"]["2.0"]),
                    *[abs(r["legs"][m] - ref["legs"][m]) for m in ("mse", "reach", "nmae", "fid", "jac")])
            ok = d < 1e-9
            ka_ok &= ok
            r["ka_vs_price246"] = dict(ref_totals=ref["totals"], ref_legs=ref["legs"], max_abs_diff=d, ok=ok)
            P(f"{'':60s}  KA vs price246 {ka_map[label]}: totals {ref['totals']['0.7']:+.5f} / {ref['totals']['1.0']:+.5f} / "
              f"{ref['totals']['2.0']:+.5f}  max|d| (totals+legs) {d:.1e} {'PASS' if ok else 'FAIL'}")
    # the difference of 246's two price rows, for comparison with the increment priced on its own ratio
    d246 = {f: P246["B2_cd4_svd10_resolved"]["totals"][f] - P246["B2_cd4_svd5_resolved"]["totals"][f] for f in ("0.7", "1.0", "2.0")}
    dlegs = {m: P246["B2_cd4_svd10_resolved"]["legs"][m] - P246["B2_cd4_svd5_resolved"]["legs"][m] for m in ("mse", "reach", "nmae", "fid", "jac")}
    P(f"{'difference of 246 price rows (svd10_res - svd5_res)':60s}{'':8s}{P246['B2_cd4_svd10_resolved']['d_pds']-P246['B2_cd4_svd5_resolved']['d_pds']:+10.6f}{'':5s}"
      f"{dlegs['mse']:+9.5f}{dlegs['reach']:+9.5f}{dlegs['nmae']:+9.5f}{dlegs['fid']:+9.5f}{dlegs['jac']:+9.5f}{'':25s}"
      f"{d246['0.7']:+9.5f}{d246['1.0']:+8.5f}{d246['2.0']:+8.5f}")
    P(f"\ngates: shape-only {GATE_SHAPE:+.4f}   amplitude-changing {GATE_AMPL:+.4f}   PLOLO pds LOO sigma {SIGMA_PDS} raw   "
      f"ordering noise sd between digest-identical arms {NOISE_ORDERING} overall")
    inc = out["INCREMENT: cd4svd10 vs cd4svd5 (the slot's one variable)"]
    a10 = out["cd4svd10 = B2_cd4_svd10_resolved vs h1union (KA)"]
    P(f"\nthe INCREMENT against the shape-only gate: {inc['totals']['0.7']/GATE_SHAPE:+.2f}x / {inc['totals']['1']/GATE_SHAPE:+.2f}x / "
      f"{inc['totals']['2']/GATE_SHAPE:+.2f}x at pds transfer 0.7 / 1.0 / 2.0; the pds leg alone {inc['pds']['0.7']:+.5f} / {inc['pds']['1']:+.5f} / {inc['pds']['2']:+.5f} "
      f"against the five cosine legs' sum {sum(inc['legs'].values()):+.5f}; leave-top-3-out of the increment "
      f"{100*S['concentration']['increment']['lt3o']:+.1f}%")
    P(f"the k=10 ARM against h1union: {a10['totals']['0.7']/GATE_SHAPE:.2f}x / {a10['totals']['1']/GATE_SHAPE:.2f}x / "
      f"{a10['totals']['2']/GATE_SHAPE:.2f}x the shape-only gate (cd4svd5: {out['cd4svd5  = B2_cd4_svd5_resolved vs h1union (KA)']['totals']['0.7']/GATE_SHAPE:.2f}x / "
      f"{out['cd4svd5  = B2_cd4_svd5_resolved vs h1union (KA)']['totals']['1']/GATE_SHAPE:.2f}x / {out['cd4svd5  = B2_cd4_svd5_resolved vs h1union (KA)']['totals']['2']/GATE_SHAPE:.2f}x); "
      f"leave-top-3-out {100*S['concentration']['arm10_vs_stack']['lt3o']:.1f}% (k=5 {100*S['concentration']['arm5_vs_stack']['lt3o']:.1f}%)")
    ex10 = six246.get("B2_cd4_svd10", {})
    ex5 = six246.get("B2_cd4_svd5_resolved", {})
    P(f"\nexact `nmae` identity (six246, a magnitude bound, not a sign): B2_cd4_svd10 (FIXED sizing; the resolved k=10 row was not run through six246) "
      f"{ex10.get('overall_nmae_vs_inc', float('nan')):+.5f} (e {ex10.get('e_vs_inc', float('nan')):.3f}, delivered ratio vs a15 {ex10.get('delivered_norm_ratio')}); "
      f"B2_cd4_svd5_resolved {ex5.get('overall_nmae_vs_inc', float('nan')):+.5f} (e {ex5.get('e_vs_inc', float('nan')):.3f}, {ex5.get('delivered_norm_ratio')})")
    P(f"\nboard raw pds implied for the INCREMENT at 0.7 / 1.0 / 2.0x: {0.7*inc['d_pds']:+.5f} / {inc['d_pds']:+.5f} / {2*inc['d_pds']:+.5f} "
      f"on top of cd4svd5's; for the k=10 ARM vs h1union's 0.83572: {0.7*a10['d_pds']:+.5f} / {a10['d_pds']:+.5f} / {2*a10['d_pds']:+.5f}")
    (WORK / "price255.json").write_text(json.dumps(dict(arms=out, diff_of_246_rows=dict(totals=d246, legs=dlegs), ka_all=bool(ka_ok),
                                                        gates=dict(shape=GATE_SHAPE, amplitude=GATE_AMPL),
                                                        noise_ordering=NOISE_ORDERING), indent=1))
    (WORK / "price255.log").write_text("\n".join(lines) + "\n")
    print(f"\nKA {'PASS' if ka_ok else 'FAIL'}; wrote {WORK/'price255.json'}")
    raise SystemExit(0 if ka_ok else 1)


if __name__ == "__main__":
    main()
