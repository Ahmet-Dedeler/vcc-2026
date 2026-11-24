"""216 — the six-member price at every candidate operating point, and the gate verdicts.

Joins the two priceable routes into one table:
  `pds` and `mse`   from `price216.py` (PLOLO x1.011; 096/197's `X`-only identity)
  `nmae`            from `six216.py`   (188's EXACT identity, wald89 pool) -- NOT 187's
                    proxy, which `price216.py` uses for the wide grid and which is ~0 by
                    construction there
  `reach` / `fid`   ordering half from `six216.py` (favourable, unpriced); yield halves
                    unreadable (rule 9), h1-test vacuous (`vacuity216.py`)
  `jac`             unpriced and NOT exempt -- the call set moves with the profile

Rule 13: this is a PROFILE-class arm family (no emission-law change), so no `reach` toll
(169 falsified it at +0.0006 for this class).  Gate +0.0019 shape-only, +0.0089 if amplitude
moves.  **211's standing rule: `claim` is at 100.0000% of its cap in the binding context, so
no arm may be re-priced on amplitude grounds** -- which is why the library arms, whose alpha
must rise 1.6-4.2x, take the amplitude gate and get no amplitude credit.

Known-answer rows, all from 207 §8.1's own selection table:
  the six-member price must read +0.00437 at f = 0.34, +0.00413 at f = 0.28 and +0.00408 at
  f = 0.40.

    uv run python final216.py
"""
from __future__ import annotations

import json
from pathlib import Path

WORK = Path(__file__).resolve().parent / "work"
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089
ALPHA0 = 0.982

SIX_KEY = {
    0.10: "a15 f=0.10 SHIPPED", 0.20: "res f=0.20", 0.28: "res f=0.28",
    0.34: "res f=0.34 <- INCUMBENT", 0.40: "res f=0.40", 0.55: "res f=0.55",
    0.80: "res f=0.80", 1.00: "res f=1.00 <- cos argmax",
}
LIB_KEY = {
    "lib w=1 (1/k) @a0": "lib1@a0",
    "lib w=1 (1/k) @adel": "lib1@adel",
    "lib w=0.5 @adel": "lib0.5@adel",
    "lib w=0.5 + res f=0.34 @adel": "lib0.5+res0.34@adel",
}
EXPECT_207 = {0.34: 0.00437, 0.28: 0.00413, 0.40: 0.00408}


def main() -> None:
    price = {r["f"]: r for r in json.loads((WORK / "price216.json").read_text())["rows"]}
    six = json.loads((WORK / "six216.json").read_text())
    form = {r["arm"]: r for r in json.loads((WORK / "formanalyse216.json").read_text())["rows"]}

    rows, gates = [], []
    print(f"{'operating point':<30}{'ovr pds':>10}{'ovr mse':>10}{'ovr nmae':>11}"
          f"{'TOTAL':>10}{'TOTAL cal':>11}{'d prec':>10}{'alpha':>8}{'class':>11}")
    print("-" * 111)
    for f, key in SIX_KEY.items():
        p, s = price[f], six[key]
        tot = p["o_pds"] + p["o_mse_k1"] + s["overall_nmae_vs_a15"]
        tot_c = p["o_pds"] + p["o_mse_kcal"] + s["overall_nmae_vs_a15"]
        rows.append(dict(name=f"res f={f:.2f}", f=f, o_pds=p["o_pds"], o_mse=p["o_mse_k1"],
                         o_mse_cal=p["o_mse_kcal"], o_nmae=s["overall_nmae_vs_a15"],
                         total_k1=tot, total_kcal=tot_c,
                         d_precision=s["d_precision_vs_a15"], alpha=ALPHA0,
                         delivered_ratio=s["delivered_norm_ratio"],
                         pds_sigma_vs_shipped=p["pds_sigma"], arm_class="shape"))
        print(f"{'res f=%.2f' % f:<30}{p['o_pds']:>+10.5f}{p['o_mse_k1']:>+10.5f}"
              f"{s['overall_nmae_vs_a15']:>+11.6f}{tot:>+10.5f}{tot_c:>+11.5f}"
              f"{s['d_precision_vs_a15']:>+10.5f}{ALPHA0:>8.3f}{'shape':>11}")
        if f in EXPECT_207:
            ok = abs(round(tot, 5) - EXPECT_207[f]) <= 6e-5
            gates.append(dict(name=f"207 §8.1 six-member price at f={f:g}",
                              expected=EXPECT_207[f], measured=round(tot, 5),
                              passed=bool(ok)))

    for key, farm in LIB_KEY.items():
        s = six[key]
        fr = form.get(farm)
        if fr is None:
            continue
        o_pds = fr["o_pds"]
        o_mse = fr["o_mse"] if fr["o_mse"] is not None else None
        amp = abs(s["alpha"] / ALPHA0 - 1.0) > 0.02
        tot = (o_pds + (o_mse or 0.0) + s["overall_nmae_vs_a15"]) if o_mse is not None else None
        rows.append(dict(name=key, f=None, o_pds=o_pds, o_mse=o_mse,
                         o_nmae=s["overall_nmae_vs_a15"], total_k1=tot,
                         d_precision=s["d_precision_vs_a15"], alpha=s["alpha"],
                         alpha_ratio=s["alpha"] / ALPHA0,
                         arm_class="AMPLITUDE" if amp else "shape"))
        t = f"{tot:>+10.5f}" if tot is not None else f"{'--':>10}"
        m = f"{o_mse:>+10.5f}" if o_mse is not None else f"{'--':>10}"
        print(f"{key:<30}{o_pds:>+10.5f}{m}{s['overall_nmae_vs_a15']:>+11.6f}{t}"
              f"{'':>11}{s['d_precision_vs_a15']:>+10.5f}{s['alpha']:>8.3f}"
              f"{('AMPLITUDE' if amp else 'shape'):>11}")

    inc = [r for r in rows if r.get("f") == 0.34][0]
    print(f"\ndeltas against the incumbent `cd4f034` (res f=0.34), which is "
          f"{inc['total_k1']:+.5f} against the shipped f=0.10:")
    print(f"{'candidate':<30}{'d overall':>12}{'gate':>9}{'x gate':>9}  verdict")
    print("-" * 76)
    verdicts = []
    for r in rows:
        if r is inc or r["total_k1"] is None:
            continue
        d = r["total_k1"] - inc["total_k1"]
        gate = GATE_AMP if r["arm_class"] == "AMPLITUDE" else GATE_SHAPE
        v = "SHIP" if d > gate else ("VETO" if d <= gate else "?")
        verdicts.append(dict(name=r["name"], delta=d, gate=gate, verdict=v,
                             arm_class=r["arm_class"]))
        print(f"{r['name']:<30}{d:>+12.5f}{gate:>9.4f}{(d/gate):>9.2f}  {v}")

    best = max(verdicts, key=lambda v: v["delta"])
    out = dict(rows=rows, gates=gates, verdicts=verdicts,
               incumbent_total_vs_shipped=inc["total_k1"],
               best_candidate=best,
               anything_clears=bool(best["delta"] > best["gate"]))
    print(f"\nbest candidate: {best['name']}  {best['delta']:+.5f} against a "
          f"{best['gate']:+.4f} gate -> {'SHIP' if out['anything_clears'] else 'VETO'}")
    print(f"known-answer gates vs 207 §8.1: "
          f"{sum(g['passed'] for g in gates)}/{len(gates)} pass  "
          f"{[(g['name'].split('=')[-1], g['measured'], g['expected']) for g in gates]}")
    (WORK / "final216.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
