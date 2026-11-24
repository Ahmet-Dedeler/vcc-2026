"""317 — is `reach` an operating-point property or a model property? Measured on the board.

094 regressed each raw member on `log(raw mse)` across the top 60 teams and found `reach`
R^2 = 0.097 with our own residual at −0.15 sd. This re-runs that on the live board and adds the
two facts that decide whether 094's "the term is already saturated" is a property of the METRIC
or of OUR CONSTRUCTION.
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SNAP = ROOT / "data/leaderboard/2026-09-07T1853Z.json"
ent = json.loads(SNAP.read_text())["entries"]
K = {"mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
     "reach": "de_wilcoxon_direction_reach_raw", "pds": "pds_cosine",
     "fid": "de_wilcoxon_direction_fidelity_yield_raw", "jac": "de_wilcoxon_sig_jaccard"}
ok = [e for e in ent if all(e.get(v) is not None for v in K.values())
      and e.get("score_avg") is not None]
ok.sort(key=lambda e: -e["score_avg"])
print(f"[317] {len(ok)} scored rows on {SNAP.name}")

print("\n[317] the `reach` scale, measured, not assumed")
r = np.array([e[K["reach"]] for e in ok])
print(f"  metric baseline b = 0.079015 (fitted), replicate anchor r = b + span = "
      f"{0.079015 + 0.885829:.6f}   (AGENTS.md §4: 0.958-0.978)")
print(f"  field raw `reach`: min {r.min():.5f}  median {np.median(r):.5f}  "
      f"p90 {np.percentile(r, 90):.5f}  max {r.max():.5f}")
us = [e for e in ok if abs(e["score_avg"] - 0.2177189) < 1e-6]
ours = us[0][K["reach"]] if us else 0.27428
print(f"  ours {ours:.5f} = {(ours - 0.079015) / 0.885829 * 100:.1f}% of the way from the "
      f"baseline to the replicate anchor; the leader {r.max():.5f} = "
      f"{(r.max() - 0.079015) / 0.885829 * 100:.1f}%")

for N in (60, 200, len(ok)):
    top = ok[:N]
    x = np.log(np.array([e[K["mse"]] for e in top]))
    print(f"\n[317] top {N}: each raw member regressed on log(raw mse)")
    for m in ("reach", "nmae", "fid", "jac", "pds"):
        y = np.array([e[K[m]] for e in top])
        sl, ic = np.polyfit(x, y, 1)
        pred = sl * x + ic
        r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        print(f"    {m:<6} slope {sl:+8.4f}/log-unit   R2 {r2:6.3f}   "
              f"residual sd {np.std(y - pred):.4f}")
