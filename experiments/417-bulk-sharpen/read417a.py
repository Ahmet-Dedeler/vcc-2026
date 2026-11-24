"""Read 417a's board row the way 420's T-read table does: implied sharpening transfer k_sh from pds and from mse.
usage: uv run python experiments/417-bulk-sharpen/read417a.py <entry id>"""
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
m = joint187.Model(joint187.F_FIT)
claim = {c: m.t[c]["claim"] for c in CTX}
Kc = {c: m.t[c]["Kc"] for c in CTX}
sumD = joint187.sumD0


def find(x):
    if isinstance(x, dict):
        if "score_avg" in x:
            return x
        for v in x.values():
            y = find(v)
            if y:
                return y


def solve_X(AB, raw):
    num = np.mean([AB[c]["A"] + Kc[c] - min(claim[c], AB[c]["B"]) for c in CTX])
    return float((num - (raw - 1.0) * sumD) / 2.0)


row = find(json.loads((ROOT / f"data/board-entries/{sys.argv[1]}.json").read_text()))
raw_pds, raw_mse = row["pds_cosine"], row["expr_mse_unbiased_capped_norm"]
dec = json.loads((ROOT / "experiments/414-bulk-reprice/work/decomp414.json").read_text())
cc = json.loads((ROOT / "experiments/414-bulk-reprice/work/cliffcheck414.json").read_text())["the-pigeon-passed-peer-review.h5ad"]
AB = {c: dict(A=cc["per_context"][c]["A"], B=cc["per_context"][c]["B"]) for c in CTX}
Xb = solve_X(dec["410b"], 0.730994)
X = solve_X(AB, raw_mse)
k_mse = (X / (Xb * 0.9798) - 1) / 0.1013
k_pds = (raw_pds - 0.851156) / 0.01713 / 0.895
print(f"417a raw pds {raw_pds:.6f} (410b 0.851156)  raw mse {raw_mse:.6f} (410b 0.730994)  overall {row['score_avg']:.7f}")
print(f"implied sharpening transfer: k_sh(pds) {k_pds:+.3f}   k_sh(mse) {k_mse:+.3f}   |diff| {abs(k_pds - k_mse):.3f} "
      f"({'CONSISTENT' if abs(k_pds - k_mse) <= 0.3 else 'DISAGREE > 0.3: the cliff model failed; re-price re-level rows'})")
print(f"X 417a {X:.6f} vs 410b {Xb:.6f} -> board relative cosine change {X / (Xb * 0.9798) - 1:+.4f} (PLOLO cancer +0.1013)")
