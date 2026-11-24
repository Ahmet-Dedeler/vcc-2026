"""800 — per-class over-read factor of the fresh (arc800) panels against an Arc-drawn panel, with a CI (coordinator's
request 10-16 10:45Z, item 3).

For each arm and held line h: v_L(k, h) on the fresh folds Lk/<h> (k = 1..6) and v_A(h) on an A/B/C-panel fold of the same
held line (`--abc ABCp`, or `LABC` once 801 registers it: A/B/C's panel on the fresh-fold recipe, the like-for-like
comparison). The class factor pools arms and held lines as a ratio of sums:
    F = Σ_arms Σ_h mean_k v_L(k, h) / Σ_arms Σ_h v_A(h)
Two-level bootstrap (4,000 draws, seed 800): resample the six panel draws with replacement (one index set shared by the
three held lines, because Lk is one panel draw across them), and resample targets with replacement inside every fold.
The panel draw dominates the spread; the target level adds the within-panel noise. A ratio whose denominator's interval
reaches 0 is reported as unbounded.

Per-fold value = stats800's fold overall (literal pds + mse via 518's identity), board units. Light: MB arrays.
    uv run python experiments/800-checker/overread800.py [--abc ABCp|LABC]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panelrank800 as P  # noqa: E402  (its roots and arm table)
import stats800 as S  # noqa: E402

X = P.X
ROOTS_LABC = [X / "801-lab-folds/work/ka-labc/bulk", X / "810-strong-direction/work/labc/bulk",
              X / "811-lolo-model/work/run-m811labc/bulk"]
CLASSES = {
    "key": ["495 cap6-cap3 (key)", "705 key rider - inc"],
    "direction": ["707 direction - inc", "sgu1 direction - inc", "811 c0_l25 - inc"],
    "shape": ["479p2g cap3+P2g-cap3 (shape)", "479sp2 cap3+S_P2-cap3 (shape)"],
}
MEMBERS = ("pds", "mse")
FRESH = "L"


def load_pair(d: Path, cand: str, base: str):
    b, c = d / f"{base}.npz", d / f"{cand}.npz"
    if not (b.exists() and c.exists()):
        return None
    return dict(np.load(b, allow_pickle=True)), dict(np.load(c, allow_pickle=True))


def fold_overall(pair, idx=None) -> float:
    zb, zc = pair
    raw = S.bulk_fold_raw(zb, zc, idx)
    return S._fold_overall(raw, "bulk", MEMBERS)


def collect(abc_tag: str):
    """{arm: {held: {"L": {k: pair}, "A": pair}}} for the arms present on both sides."""
    out = {}
    for arm, (cand, base, roots, _) in P.PAIRS.items():
        roots = list(roots) + (ROOTS_LABC if abc_tag == "LABC" else []) + (P.ROOTS_A if FRESH == "A" else [])
        per = {}
        for h in P.HELD:
            L, A = {}, None
            for r in roots:
                for d in sorted(r.glob(f"*__{h}")):
                    tag = d.name.split("__")[0]
                    p = load_pair(d, cand, base)
                    if p is None:
                        continue
                    if tag[:1] == FRESH and tag[1:].isdigit():
                        L.setdefault(int(tag[1:]), p)
                    elif tag == abc_tag and A is None:
                        A = p
            if len(L) >= 3 and A is not None:
                per[h] = dict(L=L, A=A)
        if per:
            out[arm] = per
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abc", default="ABCp", choices=["ABCp", "ABC", "LABC"])
    ap.add_argument("--boot", type=int, default=4000)
    ap.add_argument("--fresh", default="L", choices=["L", "A"], help="fresh set: L1-L6 (arc800) or A1-A6 (A-only)")
    a = ap.parse_args()
    global FRESH
    FRESH = a.fresh
    data = collect(a.abc)
    rng = np.random.default_rng(800)
    panels = sorted({k for arm in data.values() for h in arm.values() for k in h["L"]})
    res = {}
    for cls, arms in CLASSES.items():
        arms = [x for x in arms if x in data]
        if not arms:
            continue
        cells = [(arm, h) for arm in arms for h in data[arm]]
        num = sum(np.mean([fold_overall(data[arm][h]["L"][k]) for k in data[arm][h]["L"]]) for arm, h in cells)
        den = sum(fold_overall(data[arm][h]["A"]) for arm, h in cells)
        draws_n, draws_d = np.empty(a.boot), np.empty(a.boot)
        for b in range(a.boot):
            ks = rng.choice(panels, size=len(panels), replace=True)
            n_sum = d_sum = 0.0
            for arm, h in cells:
                Ld, A = data[arm][h]["L"], data[arm][h]["A"]
                vals = []
                for k in ks:
                    if k in Ld:
                        p = Ld[k]
                        vals.append(fold_overall(p, rng.integers(0, p[0]["d"].size, p[0]["d"].size)))
                n_sum += np.mean(vals) if vals else 0.0
                d_sum += fold_overall(A, rng.integers(0, A[0]["d"].size, A[0]["d"].size))
            draws_n[b], draws_d[b] = n_sum, d_sum
        ratio = draws_n / draws_d
        d_lo, d_hi = np.percentile(draws_d, [2.5, 97.5])
        bounded = d_lo > 0 or d_hi < 0
        rec = dict(arms=arms, cells=len(cells), fresh_sum=float(num), abc_sum=float(den),
                   factor=float(num / den) if den else None,
                   factor_ci95=[float(x) for x in np.percentile(ratio, [2.5, 97.5])] if bounded else None,
                   abc_sum_ci95=[float(d_lo), float(d_hi)], fresh_sum_ci95=[float(x) for x in np.percentile(draws_n, [2.5, 97.5])],
                   note=None if bounded else "A/B/C-panel sum's interval reaches 0: the factor is unbounded")
        res[cls] = rec
        ci = rec["factor_ci95"]
        print(f"{cls:9s} vs {a.abc:5s} arms {len(arms)} cells {len(cells):2d}: fresh {num:+.5f} / A/B/C-panel {den:+.5f} = "
              f"factor {rec['factor']:.2f} " + (f"[{ci[0]:.2f}, {ci[1]:.2f}]" if ci else "[unbounded]")
              + f" | A/B/C-panel sum CI [{d_lo:+.5f}, {d_hi:+.5f}]", flush=True)
    (HERE / "work").mkdir(exist_ok=True)
    (HERE / "work" / (f"overread800-{a.abc}.json" if FRESH == "L" else f"overread800-{a.abc}-fresh{FRESH}.json")).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
