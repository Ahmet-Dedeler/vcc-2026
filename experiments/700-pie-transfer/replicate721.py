#!/usr/bin/env python3
"""Independent replication of 721 §2 (scaling constants) and the DE-opportunity table, from 721's work/entries_all.json,
with code written without reading 721's scripts ([020377]'s request). Our reference row = our best entry in the file."""
import json, numpy as np
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
E = json.load(open(ROOT / "experiments/721-board-as-instrument/work/entries_all.json"))
pairs = {"fid": ("de_wilcoxon_direction_fidelity_yield_raw", "score_fid"), "jac": ("de_wilcoxon_sig_jaccard", "score_jac"),
         "reach": ("de_wilcoxon_direction_reach_raw", "score_reach"), "pds": ("pds_cosine", "score_pds"),
         "nmae": ("de_wilcoxon_lfc_nmae", "score_nmae"), "mse": ("expr_mse_unbiased_capped_norm", "score_mse")}
print("scaling constants (robust LS on entries away from clamps): score = (raw - b) / (r - b)")
for m, (rk, sk) in pairs.items():
    x = np.array([e[rk] for e in E if e.get(rk) is not None and e.get(sk) is not None], float)
    y = np.array([e[sk] for e in E if e.get(rk) is not None and e.get(sk) is not None], float)
    ok = np.isfinite(x) & np.isfinite(y) & (np.abs(y) < 3)
    if m == "mse": ok &= (y > 0.001) & (y < 0.999)          # clamp [0, 1]
    if m == "nmae": ok &= (y > -5.9)                         # floor -6
    A = np.vstack([x[ok], np.ones(ok.sum())]).T
    coef, res, *_ = np.linalg.lstsq(A, y[ok], rcond=None)
    k, c = coef; b = -c / k; r = b + 1 / k
    resid = y[ok] - A @ coef
    print(f"  {m:5s} b {b:.5f} r {r:.5f}  n {ok.sum():6d}  rms resid {np.sqrt(np.mean(resid**2)):.2e}")
# DE opportunity at equal pds: latest entry per team; teams within +-0.01 raw pds of our best entry
ours = max((e for e in E if "Dedeler" in (e.get("team") or "")), key=lambda e: e["overall"])
latest = {}
for e in E:
    t = e.get("team")
    if t and (t not in latest or e["date"] > latest[t]["date"]):
        latest[t] = e
peers = [e for t, e in latest.items() if "Dedeler" not in t and e.get("pds_cosine") is not None
         and abs(e["pds_cosine"] - ours["pds_cosine"]) <= 0.01]
print(f"our reference: {ours['model'][:40]!r} {ours['date'][:10]} overall {ours['overall']:.4f} pds {ours['pds_cosine']:.4f}; peers n={len(peers)}")
tot = 0.0
for m in ("fid", "reach", "nmae", "jac", "mse"):
    sk = pairs[m][1]
    med = float(np.median([e[sk] for e in peers]))
    d = (med - ours[sk]) / 6
    if m != "mse": tot += d
    print(f"  {m:5s} median {med:+.4f} ours {ours[sk]:+.4f} -> {d:+.4f} overall")
print(f"  DE members at the field median for our bulk: {tot:+.4f} overall")
