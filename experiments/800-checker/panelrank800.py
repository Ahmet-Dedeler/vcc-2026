"""800 — is A/B/C's panel inside the spread of the arc800 panel draws? (lead, 10-16; LAB-LOG 09:36Z open question)

For one held truth line, the fresh folds L1–L6/<held> and the A/B/C-panel folds ABC/<held> (062 operator) and ABCp/<held>
(day operator, the L folds' recipe) price the same candidate against the same truth; only the 300-target panel differs
(and, for ABC/, the operator). If A/B/C's panel were one more draw from arc800's rule, ABCp's per-fold value would be
exchangeable with the six L values: its rank among the seven is uniform on 1..7. A rank stuck at the bottom across arms
and held lines says arc800's draw rule misses something about how Arc builds panels.

Reads only arrays already on disk (standard format, bulk/<fold>/<arm>.npz). Per-fold value = stats800's fold overall
(literal pds + mse via 518's identity), board units. Light: MB arrays, seconds.
    uv run python experiments/800-checker/panelrank800.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import stats800 as S  # noqa: E402

X = ROOT / "experiments"
R804 = X / "804-known-answers/work/run/bulk"
RP = X / "801-lab-folds/work/ka-abcp/bulk"
R810, R810P = X / "810-strong-direction/work/price-final/bulk", X / "810-strong-direction/work/abcp/bulk"
R811 = [X / f"811-lolo-model/work/{r}/bulk" for r in ("run-m811s1", "run-m811s2", "run-m811s2b")]
HELD = ("k562_gwps", "orion_hct116", "orion_hek293t")
# LABC (801 follow-up 4): A/B/C's 300 targets on the fresh-fold recipe, the like-for-like comparator; read when present
# A-only folds (801 follow-up 5: arc800 minus the essential screens; DESIGN §3.1 amendments): the corrected fresh set
ROOTS_A = [X / "801-lab-folds/work/ka-aonly/bulk", X / "810-strong-direction/work/aonly/bulk",
           X / "811-lolo-model/work/run-m811aonly/bulk"]
ARC_EXTRA = []  # --arc tags (kit-built Arc-drawn panels), compared like LABC
FRESH = "L"     # which fresh set ranks the A/B/C-panel folds: "L" (L1-L6, arc800) or "A" (A1-A6, A-only); --fresh
ROOTS_LABC = [X / "801-lab-folds/work/ka-labc/bulk", X / "810-strong-direction/work/labc/bulk",
              X / "811-lolo-model/work/run-m811labc/bulk"]
# name -> (cand, base, roots, same arms on every surface?)
PAIRS = {
    "495 cap6-cap3 (key)": ("cap6", "cap3", [R804, RP], True),
    "479p2g cap3+P2g-cap3 (shape)": ("cap3+P2g", "cap3", [R804, RP], True),
    "479sp2 cap3+S_P2-cap3 (shape)": ("cap3+S_P2", "cap3", [R804, RP], True),
    "705 key rider - inc": ("705", "inc", [R804, RP], False),
    "707 direction - inc": ("707", "inc", [R804, RP], False),
    "sgu1 direction - inc": ("sgu1", "inc", [R810, R810P], False),
    "811 c0_l25 - inc": ("c0_l25", "inc", R811, False),
}


def fold_value(d: Path, cand: str, base: str):
    b, c = d / f"{base}.npz", d / f"{cand}.npz"
    if not (b.exists() and c.exists()):
        return None
    zb, zc = dict(np.load(b, allow_pickle=True)), dict(np.load(c, allow_pickle=True))
    res = S.price([S.Fold(d.name.replace("__", "/"), "bulk", zb, zc, meta={})], members=("pds", "mse"), n_boot=50, seed=800)
    return float(res["per_fold"][0]["overall"])


def is_fresh(k: str) -> bool:
    return k[:1] == FRESH and k[1:].isdigit()


def main():
    global FRESH
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--fresh", default="L", choices=["L", "A"])
    ap.add_argument("--arc", action="append", default=[], help="extra Arc-drawn comparator built by 801's kit (e.g. DEF, DEFo, "
                    "LABCk); its arrays are read from 801 work/ka-<tag>, 810 work/kit-<TAG>, 811 work/run-m811kit-<TAG>")
    args = ap.parse_args()
    FRESH = args.fresh
    for t in args.arc:
        ARC_EXTRA.append(t)
        ROOTS_LABC.extend([X / f"801-lab-folds/work/ka-{t.lower()}/bulk", X / f"810-strong-direction/work/kit-{t}/bulk",
                           X / f"811-lolo-model/work/run-m811kit-{t}/bulk"])
    out, ranks_p, ranks_062, ranks_labc = {}, [], [], []
    ranks_extra = {t: [] for t in ARC_EXTRA}
    for name, (cand, base, roots, same) in PAIRS.items():
        rec = {}
        for h in HELD:
            v = {}
            for r in list(roots) + ROOTS_LABC + (ROOTS_A if FRESH == "A" else []):
                for d in sorted(r.glob(f"*__{h}")):
                    x = fold_value(d, cand, base)
                    if x is not None:
                        v.setdefault(d.name.split("__")[0], x)
            L = np.array([v[k] for k in sorted(v) if is_fresh(k)])
            if len(L) < 3:
                continue
            row = dict(L={k: v[k] for k in sorted(v) if is_fresh(k)}, L_mean=float(L.mean()), L_sd=float(L.std(ddof=1)))
            for tag, key, acc in (("ABCp", "ABCp", ranks_p), ("ABC062", "ABC", ranks_062), ("LABC", "LABC", ranks_labc),
                                  *((t, t, ranks_extra[t]) for t in ARC_EXTRA)):
                if key in v:
                    rank = int(1 + np.sum(L < v[key]))          # 1 = below every L fold
                    row[tag] = dict(value=v[key], rank_of=f"{rank}/{len(L) + 1}",
                                    z=float((v[key] - L.mean()) / L.std(ddof=1)) if L.std(ddof=1) > 0 else None)
                    acc.append((name, h, rank, len(L) + 1, same))
            rec[h] = row
            p, o, q = row.get("ABCp", {}), row.get("ABC062", {}), row.get("LABC", {})
            print(f"{name:32s} {h:14s} L {L.mean():+.5f} sd {L.std(ddof=1):.5f} [{L.min():+.5f}, {L.max():+.5f}] "
                  f"| ABCp {p.get('value', float('nan')):+.5f} rank {p.get('rank_of', '-'):5s} z {p.get('z') or float('nan'):+.2f} "
                  f"| ABC062 {o.get('value', float('nan')):+.5f} rank {o.get('rank_of', '-'):5s}"
                  + (f" | LABC {q['value']:+.5f} rank {q['rank_of']:5s} z {q['z'] or float('nan'):+.2f}" if q else "")
                  + "".join(f" | {t} {row[t]['value']:+.5f} rank {row[t]['rank_of']:5s}" for t in ARC_EXTRA if t in row), flush=True)
        out[name] = rec
    for tag, acc in (("ABCp", ranks_p), ("ABC062", ranks_062), ("LABC", ranks_labc), *((t, ranks_extra[t]) for t in ARC_EXTRA)):
        if not acc:
            continue
        r = np.array([a[2] for a in acc]); n = np.array([a[3] for a in acc])
        u = (r - 0.5) / n                                         # mid-rank percentile; uniform under exchangeability
        bottom = int(np.sum(r == 1))
        print(f"{tag}: {len(acc)} (arm, held) cells; mean percentile {u.mean():.2f} (0.50 if exchangeable); "
              f"bottom of its line's folds in {bottom}/{len(acc)} (expected {np.sum(1 / n):.1f})")
        out[f"summary_{tag}"] = dict(cells=len(acc), mean_percentile=float(u.mean()), bottom=bottom,
                                     expected_bottom=float(np.sum(1 / n)),
                                     same_arm_cells=[a[:3] for a in acc if a[4]])
    (HERE / "work").mkdir(exist_ok=True)
    suffix = ("" if FRESH == "L" else f"-{FRESH}") + "".join(f"-{t}" for t in ARC_EXTRA)
    (HERE / "work" / f"panelrank800{suffix}.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
