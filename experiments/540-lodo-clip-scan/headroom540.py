"""540 side check (no truth beyond PLOLO's; saved arrays only): how much of the per-target amplitude headroom does the
LODO key capture, and do simple transforms of it capture more?

From 420's saved per-target arrays (x_p = <pred_p, truth_p>, a_p = |pred_p|^2 per frame, d_p = |truth_p|^2; realised 16):
cos_agg(s) = mean_frames sum_p s_p x_p / sqrt(sum_p s_p^2 a_p * sum_p d_p) for a per-target rescaling s. First-order in s
(the render's expm1 and the ballast are not linear): validated here against comp420's exact V1 -> V1_5_3 rows (the
analytic change has the right sign and ~70 % of the size).
  oracle     s_p = max(sum_f x_pf, 0) / sum_f a_pf (the per-target regression coefficient; uses truth)
  family     w = clip(w_V1^alpha * (strength / median)^beta, lo, 3^max(alpha, 1)); strength ~ sqrt(mean_f a_p) / w_V1

    uv run python experiments/540-lodo-clip-scan/headroom540.py
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
F = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")


def cosagg(x, a, d, s):
    return float(np.mean([(s * x[j]).sum() / np.sqrt((s * s * a[j]).sum() * d.sum()) for j in range(3)]))


def main():
    out = dict(oracle={}, validation={}, family=[])
    data = {}
    for f in F:
        z = np.load(ROOT / f"experiments/420-bulk-stack/work/arrays420-{f}.npz")
        d = z["d"]
        v1 = np.load(ROOT / f"experiments/418-bulk-composition/work/arrays418-{f}.npz")["V1__PP"]
        c420 = json.loads((ROOT / f"experiments/420-bulk-stack/work/comp420-{f}.json").read_text())["arms"]
        for nm in ("stack|K3+sh1.75", "stack|V1_5_3+sh1.75"):
            x, a = z[f"{nm}@R16__x"], z[f"{nm}@R16__a"]
            one = np.ones(d.size)
            orc = np.maximum(x.sum(0), 0) / np.maximum(a.sum(0), 1e-300)
            out["oracle"][f"{f}|{nm}"] = dict(cos=cosagg(x, a, d, one), oracle=cosagg(x, a, d, orc))
        x, a = z["stack|V1+sh1.75@R16__x"], z["stack|V1+sh1.75@R16__a"]
        s = np.maximum(v1, 0.5) / v1
        an = cosagg(x, a, d, s) / cosagg(x, a, d, np.ones_like(s)) - 1
        ex = c420["stack|V1_5_3+sh1.75"]["realised"]["16"]["cos_agg"] / c420["stack|V1+sh1.75"]["realised"]["16"]["cos_agg"] - 1
        out["validation"][f] = dict(analytic_r=an, exact_r=ex)
        data[f] = (x, a, d, v1, np.sqrt(a.mean(0)) / v1)
    for alpha in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0):
        for lo in (0.3, 0.5, 0.7):
            for beta in (-0.5, -0.25, 0.0, 0.25):
                r = {}
                for f, (x, a, d, v1, st) in data.items():
                    wn = np.clip(v1 ** alpha * (st / np.median(st)) ** beta, lo, 3.0 ** max(alpha, 1))
                    r[f] = cosagg(x, a, d, wn / v1) / cosagg(x, a, d, np.maximum(v1, 0.5) / v1) - 1
                out["family"].append(dict(alpha=alpha, lo=lo, beta=beta, r=r, cancer=float(np.mean([r[f] for f in F[:3]]))))
    for k, v in out["oracle"].items():
        print(f"{k:40s} cos {v['cos']:.5f} oracle {v['oracle']:.5f} ({100 * (v['oracle'] / v['cos'] - 1):+.1f} %)")
    for f, v in out["validation"].items():
        print(f"validation {f:14s} V1 -> V1_5_3: analytic {100 * v['analytic_r']:+.3f} % exact {100 * v['exact_r']:+.3f} %")
    best = sorted(out["family"], key=lambda t: -t["cancer"])[:5]
    for b in best:
        print(f"family alpha {b['alpha']:g} lo {b['lo']:g} beta {b['beta']:+g}: cancer {100 * b['cancer']:+.2f} % "
              + " ".join(f"{f[:6]} {100 * b['r'][f]:+.2f} %" for f in F))
    (Path(__file__).resolve().parent / "work" / "headroom540.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
