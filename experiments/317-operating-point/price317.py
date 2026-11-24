"""317 — price a raw six-member delta measured on h1-test in BOARD overall points.

Two conversions, kept apart on purpose.

(1) THE BOARD'S OWN SCALING MAP. Every leaderboard row publishes both the raw member and its
    `score_*`, so `scaled = (raw - b)/(r - b)` is an exact per-member linear fit wherever the
    per-context clamps are inactive. Fitted on the 2026-09-07T1853Z snapshot (778 entries),
    restricted to each member's UNCLAMPED core, it reproduces the board's own scaled columns
    to R^2 >= 0.9978 (four members >= 0.9998). `score_avg` is the plain mean of the six
    (max |diff| 3.3e-16 over 778 rows), so d(overall) = sum_m d(scaled_m)/6.

(2) THE ARM-CLASS TRANSFER FACTOR. A raw delta measured on h1-test is NOT a board raw delta.
    `docs/BOARD-LEDGER.md` §4 and 164 §6 put the h1-test -> board `reach` factor for the
    SPREAD-POLICY arm class at 0.39-0.4277 (h1-test OVERSTATES by ~2.3-2.6x); 221 measured
    0.769x for the designed-ORDER class and found no factor at all for two others. This
    script reports factor 1.0 (the optimistic reading) and 0.4277 side by side and never
    hides which one is being quoted.
"""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
MEMBERS = ["pds", "mse", "nmae", "fid", "reach", "jac"]
RAWKEY = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm",
          "nmae": "de_wilcoxon_lfc_nmae",
          "fid": "de_wilcoxon_direction_fidelity_yield_raw",
          "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard"}


def board_slopes(snapshot=None):
    """Per-member d(scaled)/d(raw), fitted on the board's unclamped core."""
    snap = Path(snapshot) if snapshot else sorted(
        (HERE.parents[1] / "data/leaderboard").glob("2026-09-0*Z.json"))[-1]
    ent = json.loads(snap.read_text())["entries"]
    out = {}
    for m in MEMBERS:
        x = np.array([e[RAWKEY[m]] for e in ent if e.get(RAWKEY[m]) is not None
                      and e.get(f"score_{m}") is not None], float)
        y = np.array([e[f"score_{m}"] for e in ent if e.get(RAWKEY[m]) is not None
                      and e.get(f"score_{m}") is not None], float)
        sel = np.ones(x.size, bool)
        if m == "mse":                       # 690/778 sit exactly on the 0 clamp
            sel = (y > 0.02) & (y < 0.95)
            lo, hi = np.percentile(x[sel], [10, 90])
            sel &= (x >= lo) & (x <= hi)
        if m == "nmae":                      # floored at -6
            sel = (y > -1.0) & (y < 1.5)
        sl, ic = np.polyfit(x[sel], y[sel], 1)
        r2 = 1 - ((y[sel] - (sl * x[sel] + ic)) ** 2).sum() / (
            (y[sel] - y[sel].mean()) ** 2).sum()
        out[m] = {"slope": float(sl), "b": float(-ic / sl), "span": float(1 / sl),
                  "r2": float(r2), "n": int(sel.sum()),
                  "raw_lo": float(x[sel].min()), "raw_hi": float(x[sel].max())}
    return snap.name, out


def price(delta_raw, slopes, factor=1.0):
    """A raw h1-test delta -> board overall points, at one arm-class transfer factor."""
    per = {m: delta_raw.get(m, 0.0) * factor * slopes[m]["slope"] for m in MEMBERS}
    return per, sum(per.values()) / 6.0


if __name__ == "__main__":
    name, sl = board_slopes()
    print(f"[317] board scaling map, snapshot {name}")
    for m in MEMBERS:
        s = sl[m]
        print(f"  {m:<6} d(scaled)/d(raw) {s['slope']:+.6f}  b {s['b']:+.6f}  "
              f"span {s['span']:+.6f}  R2 {s['r2']:.6f}  n {s['n']:4d}  "
              f"fitted on raw in [{s['raw_lo']:.4f}, {s['raw_hi']:.4f}]")
    (HERE / "work").mkdir(exist_ok=True)
    (HERE / "work/board_slopes.json").write_text(json.dumps({"snapshot": name, "slopes": sl},
                                                            indent=1))
    print("wrote work/board_slopes.json")
