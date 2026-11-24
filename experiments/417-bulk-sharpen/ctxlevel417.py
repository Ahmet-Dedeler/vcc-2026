"""417: per-context bulk re-levelling for a SCORED ballast artifact (414 section 6, generalised).

Reads the artifact's render-space A_c/B_c from 414's cliffcheck414.json (run cliffcheck414.py on it first), solves X from
its BOARD raw mse with 187's one-equation model (claim_c/Kc_c/sumD from joint187's fitted model, exactly as analyse414's
solve_X), and prints for each context: rho_c, the cliff bulk scale sqrt(claim_c/B_c), the EV-optimal scale under 187's
cliff sd, the predicted change in raw mse, and the VCC_BULK_CTX_SCALE string. Known answer: on 410b's A_c/B_c and board
raw it must reproduce analyse414's X for 410b (0.035723).
usage: uv run python experiments/417-bulk-sharpen/ctxlevel417.py <artifact.h5ad name in cliffcheck414.json> <board raw mse>
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/187-joint-amplitude"))
sys.path.insert(0, str(ROOT / "experiments/096-mse-floor"))
sys.path.insert(0, str(ROOT / "src"))
import joint187  # noqa: E402

CTX = ("A", "B", "C")
SD187 = {"A": 0.0025, "B": 0.0004, "C": 0.0013}
MSE_SC = -1.051
m = joint187.Model(joint187.F_FIT)
claim = {c: m.t[c]["claim"] for c in CTX}
Kc = {c: m.t[c]["Kc"] for c in CTX}
sumD = joint187.sumD0


def solve_X(AB, raw):
    num = np.mean([AB[c]["A"] + Kc[c] - min(claim[c], AB[c]["B"]) for c in CTX])
    return float((num - (raw - 1.0) * sumD) / 2.0)


def f0(c, AB, s, X):
    return 1.0 + (s * s * AB[c]["A"] + Kc[c] - min(claim[c], s * s * AB[c]["B"]) - 2.0 * s * X) / sumD


dec = json.loads((ROOT / "experiments/414-bulk-reprice/work/decomp414.json").read_text())
ka = solve_X(dec["410b"], 0.730994)
assert abs(ka - 0.035723) < 5e-6, ka
cc = json.loads((ROOT / "experiments/414-bulk-reprice/work/cliffcheck414.json").read_text())[sys.argv[1]]
AB = {c: dict(A=cc["per_context"][c]["A"], B=cc["per_context"][c]["B"]) for c in CTX}
X = solve_X(AB, float(sys.argv[2]))
print(f"KA: 410b X {ka:.6f} (analyse414 0.035723) PASS | {sys.argv[1]}: X {X:.6f} at board raw {float(sys.argv[2]):.6f}")
rng = np.random.default_rng(417)
grid = np.linspace(0.90, 1.12, 2201)
scales, tot, tot_ev = {}, 0.0, 0.0
for c in CTX:
    cl = float(np.sqrt(claim[c] / AB[c]["B"]))
    dv = np.array([f0(c, AB, s, X) - f0(c, AB, 1.0, X) for s in grid])
    draws = cl * np.exp(rng.normal(0, SD187[c], 400))

    def ev(s):
        v = 0.0
        for d in draws:
            cm = d * d * AB[c]["B"]
            v += (s * s * AB[c]["A"] - min(cm, s * s * AB[c]["B"]) - 2 * s * X) - (AB[c]["A"] - min(cm, AB[c]["B"]) - 2 * X)
        return v / len(draws) / sumD
    g2 = grid[::4]
    evs = np.array([ev(s) for s in g2])
    scales[c] = float(g2[evs.argmin()])
    tot += float(dv.min()) / 3
    tot_ev += float(evs.min()) / 3
    print(f"  {c}: rho {AB[c]['B'] / claim[c]:.4f}  cliff scale {cl:.4f}  argmin {grid[dv.argmin()]:.4f} (d raw {dv.min():+.5f})  "
          f"EV-optimal {scales[c]:.4f} (d raw {evs.min():+.5f})")
print(f"to the cliffs: d raw {tot:+.5f} -> overall {MSE_SC * tot / 6:+.5f};  EV-optimal: d raw {tot_ev:+.5f} -> overall {MSE_SC * tot_ev / 6:+.5f}")
print("VCC_BULK_CTX_SCALE=" + ",".join(f"{c}={scales[c]:.4f}" for c in CTX))
