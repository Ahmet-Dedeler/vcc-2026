"""243 — ARM6's own identity terms, precomputed from its BUILT package so that the board raw
`mse` that lands on 2026-09-05 converts directly into a MEASURED cosine (Phase 2's reading,
not the floor) and a realised PLOLO->board factor, with no further build or stream.

Reads submissions/229-stack6-g035.h5ad totals streamed by totals243.py (work/totals-ARM6.npz)
and writes work/arm6_243.json plus a reading table.

    uv run python arm6_243.py
"""
from __future__ import annotations

import json

import numpy as np

from lib243 import COEFF, CTX, ROOT, SUMD, WORK, floor_cos
from invert243 import bctrl, invert, terms

R = json.loads((WORK / "rule243.json").read_text())
I = json.loads((WORK / "invert243.json").read_text())


def main() -> None:
    genes, b = bctrl()
    z = np.load(WORK / "totals-ARM6.npz", allow_pickle=True)
    assert np.array_equal(np.asarray(z["genes"]).astype(str), genes)
    meta = json.loads((ROOT / "submissions/229-stack6-g035.h5ad.meta.json").read_text())
    assert int(z["nnz"]) == meta["output"]["nnz"], "streamed nnz != meta nnz"
    per = terms(np.asarray(z["totals"]), np.asarray(z["blocks"]).astype(str), genes, b)
    inc = I["arms"]["h1union"]
    print("ARM6 (229-stack6-g035) vs h1union, delivered terms:")
    for c in CTX:
        print(f"  {c}: A {per[c]['A']:.6f} (h1union {inc['per'][c]['A']:.6f}, "
              f"x{per[c]['A']/inc['per'][c]['A']:.5f})  B {per[c]['B']:.6f}  cap {per[c]['cap']:.6f}  "
              f"leg {'claim' if per[c]['cap'] < per[c]['B'] else 'B'}")
    amp = float(np.mean([np.sqrt(per[c]["A"] / inc["per"][c]["A"]) for c in CTX]))
    print(f"  delivered amplitude vs h1union: {100*(amp-1):+.3f}%")

    cos_inc = inc["inv"]["cos"]
    rel_plolo = R["arm6"]["rel_agg_3c"]
    raw_inc = inc["board"]["raw"]["mse"]
    msa = float(np.mean([np.sqrt(per[c]["A"]) for c in CTX]))
    num6 = float(np.mean([per[c]["A"] + per[c]["Kc"] - min(per[c]["cap"], per[c]["B"]) for c in CTX]))

    def raw_for_k(k):
        """board raw mse ARM6 would post if its measured cosine were h1union's x (1 + k rel)."""
        X = cos_inc * (1 + k * rel_plolo) * np.sqrt(SUMD) * msa
        return 1.0 + (num6 - 2 * X) / SUMD

    # the price, redone with ARM6's OWN amplitude terms (rule243's ARM6 rows assumed
    # h1union's A and B; the -1.24% delivered amplitude is charged on the budget leg here)
    SPAN_MSE, SPAN_PDS = 0.950162566, 0.452
    d_pds_plolo = R["arm6"]["d_pds_3c"]
    raw_at_k0 = raw_for_k(0.0)
    print(f"\namplitude charge: at ZERO cosine change ARM6 would post raw mse {raw_at_k0:.5f} vs "
          f"h1union's {raw_inc:.5f} -> {raw_at_k0-raw_inc:+.5f} raw = "
          f"{-(raw_at_k0-raw_inc)/SPAN_MSE/6:+.5f} overall, from its -{100*(1-amp):.2f}% amplitude alone")
    print(f"\nARM6 re-priced with its own terms (the legs: mse from the identity at ARM6's A/B; "
          f"pds = k_pds x PLOLO d pds {d_pds_plolo:+.6f} / 0.452 / 6; DE-four = rate x d cos):")
    print(f"  {'model':46s}{'k_cos':>6s}{'k_pds':>6s}{'DE4':>5s}{'d cos':>10s}{'raw mse':>9s}"
          f"{'o_mse':>9s}{'o_pds':>9s}{'o_DE4':>9s}{'TOTAL':>9s}")
    models = [
        ("new-source class, point (h1union's own rates)", 0.940, 0.705, 0.239),
        ("new-source class, a15-tag k", 1.005, 0.705, 0.239),
        ("new-source class, LOW band", 0.940, 0.705, 0.173),
        ("new-source class, HIGH band", 1.005, 1.146, 0.423),
        ("coverage-linear model (cov 0.587)", R["classes"]["a23"]["cov_intercept"]
         + R["classes"]["a23"]["cov_slope"] * (176 / 300), 0.705, 0.239),
        ("dose class (if the union point was the outlier)", R["classes"]["a23"]["k_dose"], 0.705, 0.423),
    ]
    preds, price_rows = {}, []
    for label, kc, kp, de in models:
        dcos = kc * rel_plolo * cos_inc
        raw = raw_for_k(kc)
        o_mse = (raw_inc - raw) / SPAN_MSE / 6
        o_pds = kp * d_pds_plolo / SPAN_PDS / 6
        o_de4 = de * dcos
        tot = o_mse + o_pds + o_de4
        preds[label] = raw
        price_rows.append(dict(label=label, k_cos=kc, k_pds=kp, de4=de, dcos=dcos, raw=raw,
                               o_mse=o_mse, o_pds=o_pds, o_de4=o_de4, total=tot))
        print(f"  {label:46s}{kc:>6.2f}{kp:>6.2f}{de:>5.2f}{dcos:>+10.6f}{raw:>9.5f}"
              f"{o_mse:>+9.5f}{o_pds:>+9.5f}{o_de4:>+9.5f}{tot:>+9.5f}")
    print(f"  (229's six-member measured, nmae identity leg -0.004783: -0.00219; 237's bins at +-0.0015)")
    # what the FLOOR would read at the class prediction, i.e. tomorrow's amplitude bias
    raw_cls = preds["new-source class, point (h1union's own rates)"]
    k_floor = (floor_cos(raw_cls) - inc["floor"]) / (rel_plolo * inc["floor"])
    print(f"  floor-read bias: at the class prediction (raw {raw_cls:.5f}) the FLOOR reads k = "
          f"{k_floor:.3f} against a true 0.94 -- the -{100*(1-amp):.2f}% amplitude costs "
          f"{0.94-k_floor:.2f} units of k in the floor reading.  Read the MEASURED column.")

    print(f"\nreading table for tomorrow: board raw mse -> measured cos -> realised k "
          f"(= d cos / ({100*rel_plolo:.3f}% x {cos_inc:.6f}))")
    print(f"  {'raw mse':>10s}{'floor cos':>11s}{'meas cos':>10s}{'d cos':>10s}{'rel %':>8s}"
          f"{'k realised':>12s}  reading")
    grid = sorted(set(list(np.round(np.arange(0.800, 0.8226, 0.0025), 4)) + [round(v, 5) for v in preds.values()]))
    rows = []
    for raw in grid:
        r = invert(per, float(raw))
        dcos = r["cos"] - cos_inc
        k = dcos / (rel_plolo * cos_inc)
        if k < 0.75:
            note = "below every model: the union-class factor itself was optimistic"
        elif k < 1.08:
            note = "NEW-SOURCE class (0.94-1.00)"
        elif k < 1.25:
            note = "coverage-linear (1.15)"
        elif k < 1.55:
            note = "DOSE class (1.36): the union point was the outlier"
        else:
            note = "above every model"
        tag = "  <- " + ", ".join(l for l, v in preds.items() if abs(v - raw) < 5e-6) if any(
            abs(v - raw) < 5e-6 for v in preds.values()) else ""
        rows.append(dict(raw=float(raw), floor=floor_cos(float(raw)), cos=r["cos"], dcos=dcos, k=k, note=note))
        print(f"  {raw:>10.5f}{floor_cos(float(raw)):>11.6f}{r['cos']:>10.6f}{dcos:>+10.6f}"
              f"{100*dcos/cos_inc:>+8.2f}{k:>12.3f}  {note}{tag}")
    th = {k: raw_for_k(k) for k in (0.75, 1.08, 1.25, 1.55)}
    print("\n  boundaries (board raw mse): " + "  ".join(f"k={k}: {v:.5f}" for k, v in th.items()))
    out = dict(per=per, amp_vs_h1union=amp, cos_inc=cos_inc, rel_plolo_3c=rel_plolo, rows=rows,
               thresholds=th, predicted_raw=preds, price=price_rows, raw_at_k0=raw_at_k0,
               floor_k_bias=0.94 - k_floor)
    (WORK / "arm6_243.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"wrote {WORK/'arm6_243.json'}")


if __name__ == "__main__":
    main()
