"""438: is a same-axis 2-donor ("held-out Orion") operator an acceptable SURROGATE for the full panel's 3-donor operator?

efr has no 3-donor operator (a fit peaks ~13 GB, over this mandate's 10 GB cap). Its axis (union of the donors' genes
minus the panel) is the same as either held-out-Orion operator's, because both Orion libraries carry the same genes.
Here, on the three panels that HAVE the real operator (abc, real, efp), full438 was run with each surrogate as well:
  * the anchored render energies (A_c, B_c) per design, surrogate vs real operator;
  * the energy ratios the unlevelled price reads (arm / 417a analogue);
  * the resulting prices on the panel's own folds (unlevelled and levelled), surrogate energies vs real ones.
The surrogate error on these three panels is the error bar carried by efr's energy-dependent numbers.

    uv run python surr438.py          -> work/surr438.json (printed table)
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402

import numpy as np  # noqa: E402

from price438 import REF, SLOTS, WORK, Panel  # noqa: E402

ROWS = ("h2c1|K3+sh1.25", "h1c1.5|K3+sh1.25", "h2c1.5|K3+sh1.25", "h2c1.5|K3+sh1.5", "h1c1|K3+sh1.5",
        "h2c1.5|V1_5_3+sh1.25", "h3c1.5|K3+sh1.25", "h2c2|K3+sh1.25", "h0c1.5|K3+sh1.25")


def main():
    out = {}
    for panel in ("abc", "real", "efp"):
        real = json.loads((WORK / panel / "full438.json").read_text())["designs"]
        P0 = Panel(panel)
        for s in ("s1", "s2", "fold"):
            tag = f"{panel}_{s}"
            if s == "fold":
                from price438 import fold_energies
                sur = fold_energies(panel, "held", P0.claims)
                PS = Panel(panel, energy="fold")
            else:
                p = WORK / tag / "full438.json"
                if not p.exists():
                    continue
                sur = json.loads(p.read_text())["designs"]
                PS = Panel(panel, full_tag=tag)
            common = [n for n in real if n in sur]
            eB = [sur[n]["AB"][c]["B"] / real[n]["AB"][c]["B"] - 1 for n in common for c in SLOTS]
            er = [(sur[n]["AB"][c]["B"] / sur[REF]["AB"][c]["B"]) / (real[n]["AB"][c]["B"] / real[REF]["AB"][c]["B"]) - 1
                  for n in common for c in SLOTS]
            dp = {}
            for n in ROWS:
                for lv in (False, True):
                    a = np.mean([P0.price(f, n, 1.0, REF, lv)["ov"] for f in P0.folds])
                    b = np.mean([PS.price(f, n, 1.0, REF, lv)["ov"] for f in PS.folds])
                    dp[f"{n}|{'lev' if lv else 'unlev'}"] = dict(real=float(a), surrogate=float(b), d=float(b - a))
            out[tag] = dict(B_abs_err=dict(mean=float(np.mean(eB)), max_abs=float(np.max(np.abs(eB)))),
                            ratio_err=dict(mean=float(np.mean(er)), max_abs=float(np.max(np.abs(er)))), price=dp,
                            max_abs_price_d_unlev=max(abs(v["d"]) for k, v in dp.items() if k.endswith("unlev")),
                            max_abs_price_d_lev=max(abs(v["d"]) for k, v in dp.items() if k.endswith("|lev")))
            print(f"{tag}: anchored B err mean {100 * out[tag]['B_abs_err']['mean']:+.2f}% (max |.| "
                  f"{100 * out[tag]['B_abs_err']['max_abs']:.2f}%); arm/417a ratio err mean {100 * out[tag]['ratio_err']['mean']:+.3f}% "
                  f"(max {100 * out[tag]['ratio_err']['max_abs']:.3f}%); price |d| max unlevelled "
                  f"{out[tag]['max_abs_price_d_unlev']:.5f}, levelled {out[tag]['max_abs_price_d_lev']:.5f}")
            for k, v in dp.items():
                print(f"    {k:28s} real {v['real']:+.5f} surrogate {v['surrogate']:+.5f} d {v['d']:+.5f}")
    (WORK / "surr438.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
