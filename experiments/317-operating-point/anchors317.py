"""317 — the board's own per-member scaling map, refit from the live snapshots.

094 fitted `scaled = (raw - b)/(r - b)` on 780 board entries and reproduced four members to
R^2 > 0.9999. The snapshots now carry BOTH the raw member and its `score_*`, so the map is
an exact per-member linear fit and needs no assumption. It is what converts a raw six-member
delta measured on h1-test into board OVERALL points:

    d(overall) = sum_members  d(raw_m) * slope_m / 6      (score_avg is the plain mean of six)

KNOWN-ANSWER: the fitted map must reproduce every entry's own published `score_*` and their
mean must reproduce `score_avg`.
"""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SNAP = sorted((ROOT / "data/leaderboard").glob("2026-09-0*Z.json"))[-1]
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm",
       "nmae": "de_wilcoxon_lfc_nmae", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
       "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard"}
ent = json.loads(SNAP.read_text())["entries"]
print(f"[317] snapshot {SNAP.name}: {len(ent)} entries")

fits = {}
for m, rawkey in RAW.items():
    x = np.array([e[rawkey] for e in ent if e.get(rawkey) is not None
                  and e.get(f"score_{m}") is not None], dtype=float)
    y = np.array([e[f"score_{m}"] for e in ent if e.get(rawkey) is not None
                  and e.get(f"score_{m}") is not None], dtype=float)
    # `mse` is clamped to [0,1] and `nmae` floored at -6: fit only on the unclamped interior
    keep = np.ones(x.size, bool)
    if m == "mse":
        keep = (y > 1e-9) & (y < 1 - 1e-9)
    if m == "nmae":
        keep = y > -6 + 1e-9
    A = np.vstack([x[keep], np.ones(keep.sum())]).T
    slope, icpt = np.linalg.lstsq(A, y[keep], rcond=None)[0]
    pred = A @ np.array([slope, icpt])
    ss = 1.0 - ((y[keep] - pred) ** 2).sum() / ((y[keep] - y[keep].mean()) ** 2).sum()
    fits[m] = {"slope": float(slope), "intercept": float(icpt), "r2": float(ss),
               "n": int(keep.sum()), "b": float(-icpt / slope),
               "span_r_minus_b": float(1.0 / slope),
               "max_abs_resid": float(np.abs(y[keep] - pred).max())}
    print(f"  {m:<6} slope {slope:+.6f}  b {-icpt/slope:+.6f}  (r-b) {1/slope:+.6f}  "
          f"R2 {ss:.10f}  n {keep.sum():4d}  max|resid| {np.abs(y[keep]-pred).max():.2e}")

# KNOWN-ANSWER: score_avg is the plain mean of the six published score_* columns
d = []
for e in ent:
    if all(e.get(f"score_{m}") is not None for m in RAW) and e.get("score_avg") is not None:
        d.append(abs(np.mean([e[f"score_{m}"] for m in RAW]) - e["score_avg"]))
print(f"\n[KA] score_avg == mean(six score_*) on {len(d)} entries, max |diff| {max(d):.3e}")

print("\n[317] raw-delta -> OVERALL-point conversion (slope/6), i.e. what one unit of a raw "
      "member is worth on the board's overall:")
for m in RAW:
    print(f"  1.0 raw {m:<6} = {fits[m]['slope']/6:+.6f} overall")
Path(__file__).with_name("work").mkdir(exist_ok=True)
Path(__file__).with_name("work").joinpath("anchors317.json").write_text(
    json.dumps({"snapshot": SNAP.name, "fits": fits,
                "score_avg_max_diff": max(d)}, indent=1))
print("\nwrote work/anchors317.json")
