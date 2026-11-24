"""317 — where is h1-test's own amplitude optimum, and does it sit at the shipped point?

For a predicted effect delta_hat = s * d against a truth delta, the un-normalised squared error
is quadratic in s:  ||s d - delta||^2 = s^2 ||d||^2 - 2 s <d,delta> + ||delta||^2.
`expr_mse_unbiased_capped_norm` adds a jackknife correction and a cap, but over a modest range
of s the raw member is still well approximated by a parabola. Six points (base at s = 1 plus
five amplitude arms) over-determine it, so the vertex is a measurement with a residual, not an
assumption. If the vertex sits at s ~ 1, h1-test independently confirms 187's board-side
claim that the shipped amplitude is the `mse` optimum.
"""
import json
from pathlib import Path
import numpy as np
from read_agg import read

HERE = Path(__file__).resolve().parent
CO = HERE / "cloud-output"
BASE = read(HERE.parents[1] /
            "experiments/221-cloud-scorer/cloud-output/ladder8-20260904/rc_callset/run/agg_results.csv")
S = {"amp070": 0.70, "amp085": 0.85, "base": 1.00, "amp120": 1.20,
     "amp145": 1.45, "amp200": 2.00}

pts = {}
for name, s in S.items():
    if name == "base":
        pts[s] = BASE
        continue
    p = CO / f"s317-{name}/{name}/run/agg_results.csv"
    if p.exists():
        pts[s] = read(p)
if len(pts) < 3:
    print(f"[317] only {len(pts)} amplitude points collected; need >= 3"); raise SystemExit

xs = np.array(sorted(pts))
print(f"[317] {len(xs)} amplitude points on h1-test: s = {list(xs)}")
for m in ("mse", "nmae", "reach", "fid", "jac", "pds"):
    ys = np.array([pts[s][m] for s in xs])
    print(f"\n  {m}: " + "  ".join(f"s={s:.2f} {y:.6f}" for s, y in zip(xs, ys)))
    if len(xs) >= 3:
        a, b, c = np.polyfit(xs, ys, 2)
        resid = ys - (a * xs ** 2 + b * xs + c)
        if abs(a) > 1e-12:
            vtx = -b / (2 * a)
            kind = "minimum" if a > 0 else "maximum"
            print(f"    parabola vertex ({kind}) at s = {vtx:.4f};  "
                  f"max |residual| {np.abs(resid).max():.2e}")
        # first difference around the shipped point
    if 1.00 in pts and 1.20 in pts and 0.85 in pts:
        dn = (pts[1.00][m] - pts[0.85][m]) / 0.15
        up = (pts[1.20][m] - pts[1.00][m]) / 0.20
        print(f"    slope below s=1: {dn:+.6f} per unit s;  above: {up:+.6f}")
