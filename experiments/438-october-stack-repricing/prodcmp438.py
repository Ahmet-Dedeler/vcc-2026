"""438: the day-of shortcut -- PLOLO folds on the panel's 3-donor operator ('prod', leaky: the held library is in its
basis) instead of three held-out fits ('held'). Same folds, same designs, same level factors; only the fold operator
differs. Reports the price difference per row and whether the day-of rule's choice changes.

    uv run python prodcmp438.py [panel ...]          -> work/prodcmp438.{json,md}
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402

import numpy as np  # noqa: E402

from price438 import REF, WORK, Panel  # noqa: E402

ROWS = [("H1 x2 half", "h2c1|K3+sh1.25", REF), ("CD4 x1.5 half", "h1c1.5|K3+sh1.25", REF),
        ("stack", "h2c1.5|K3+sh1.25", REF), ("CD4 x1.5 on top of H1 x2", "h2c1.5|K3+sh1.25", "h2c1|K3+sh1.25"),
        ("H1 x2 -> x3 (CD4 x1)", "h3c1|K3+sh1.25", "h2c1|K3+sh1.25"), ("q 1.25 -> 1.5 on the stack", "h2c1.5|K3+sh1.5", "h2c1.5|K3+sh1.25"),
        ("key K3 -> LODO (0.5, 3) on the stack", "h2c1.5|V1_5_3+sh1.25", "h2c1.5|K3+sh1.25")]


def main(panels):
    out, md = {}, ["| panel | step | held-out fold operators (lev k=1) | 3-donor operator on the folds | d | n+ held / prod |",
                   "|---|---|---|---|---:|---|"]
    for t in panels:
        H, Pp = Panel(t, "held"), Panel(t, "prod")
        if len(Pp.folds) < 3:
            print(f"skip {t}: prod folds {Pp.folds}")
            continue
        for lab, a, r in ROWS:
            vh = [H.price(f, a, 1.0, r, True)["ov"] for f in H.folds]
            vp = [Pp.price(f, a, 1.0, r, True)["ov"] for f in Pp.folds]
            out.setdefault(t, {})[lab] = dict(held=vh, prod=vp)
            md.append(f"| {t} | {lab} | {np.mean(vh):+.5f} | {np.mean(vp):+.5f} | {np.mean(vp) - np.mean(vh):+.5f} | "
                      f"{sum(v > 0 for v in vh)}/3 / {sum(v > 0 for v in vp)}/3 |")
        print("\n".join(md[-len(ROWS):]), flush=True)
        sys.path.insert(0, str(WORK.parent))
        import recommend438
        for ops in ("held", "prod"):
            r = recommend438.main(t, ops)
            out[t][f"choice_{ops}"] = r["choice"]
    (WORK / "prodcmp438.json").write_text(json.dumps(out, indent=1, default=lambda o: list(o) if isinstance(o, tuple) else str(o)))
    md.append("")
    for t in out:
        md.append(f"- {t}: rule choice held {out[t]['choice_held']} / prod {out[t]['choice_prod']}")
    (WORK / "prodcmp438.md").write_text("\n".join(md) + "\n")
    print("\n".join(md[-len(out):]))


if __name__ == "__main__":
    main(sys.argv[1:] or ["efp", "abc"])
