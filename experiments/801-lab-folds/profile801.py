"""801 -- truth-strength profile per fold (follow-up 4): from 804's standard-format arrays, the cap3 arm's `d` (per-target
truth energy on the held line's scored genes), per fold: median d, mean d, the top-30 share of sum d, the count of targets
with d > 4x the fold median; and the 495 per-fold gain (ka804's `cap6` arm vs its base cap3: r = relative aggregate cosine,
dpds literal). Spearman across folds between each statistic and the 495 r: pooled (all folds) and within each held line.
Folds: L1-L6 and ABC062 (= 804's ABC/<held>) from experiments/804-known-answers/work/run; ABCp from work/ka-abcp; LABC from
work/ka-labc.          uv run python experiments/801-lab-folds/profile801.py      -> work/profile801.{json,txt}
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN804 = ROOT / "experiments/804-known-answers/work/run"
HELD = ("k562_gwps", "orion_hct116", "orion_hek293t")


def fold_dirs():
    for h in HELD:
        for k in range(1, 7):
            yield f"L{k}/{h}", RUN804, f"L{k}__{h}"
        yield f"ABC062/{h}", RUN804, f"ABC__{h}"
        yield f"ABCp/{h}", HERE / "work/ka-abcp", f"ABCp__{h}"
        yield f"LABC/{h}", HERE / "work/ka-labc", f"LABC__{h}"


def main():
    rows = []
    for fid, root, stem in fold_dirs():
        z = np.load(root / "bulk" / stem / "cap3.npz")
        d = np.asarray(z["d"], float)
        m = json.loads((root / "meta" / f"{stem}.json").read_text())
        vb = m["arms"]["cap6"]["vs_base"]
        assert vb["base"] == "cap3", fid
        med = float(np.median(d))
        top30 = float(np.sort(d)[::-1][:30].sum() / d.sum())
        rows.append(dict(fold=fid, held=fid.split("/")[1], n=int(d.size), median_d=med, mean_d=float(d.mean()),
                         top30_share=top30, n_gt4med=int((d > 4 * med).sum()), mean_over_median=float(d.mean() / med),
                         r495=float(vb["r"]), dpds495=float(vb["dpds"])))
    stats = ("median_d", "mean_d", "top30_share", "n_gt4med", "mean_over_median")
    corr = {}
    for scope in ("pooled",) + HELD:
        R = rows if scope == "pooled" else [r for r in rows if r["held"] == scope]
        corr[scope] = {}
        for s in stats:
            for g in ("r495", "dpds495"):
                rho, p = spearmanr([r[s] for r in R], [r[g] for r in R])
                corr[scope][f"{s}~{g}"] = dict(rho=float(rho), p=float(p), n=len(R))
    (HERE / "work/profile801.json").write_text(json.dumps(dict(rows=rows, spearman=corr), indent=1))
    out = ["fold               n  median_d   mean_d  top30  n>4med  mean/med   495 r     495 dpds"]
    for r in rows:
        out.append(f"{r['fold']:18s} {r['n']:3d} {r['median_d']:9.4f} {r['mean_d']:8.4f} {r['top30_share']:6.3f} {r['n_gt4med']:6d} "
                   f"{r['mean_over_median']:8.3f} {100 * r['r495']:+8.3f}% {r['dpds495']:+9.5f}")
    out.append("")
    for scope, c in corr.items():
        out.append(f"Spearman ({scope}, n={c['median_d~r495']['n']}): " + "; ".join(
            f"{k} {v['rho']:+.2f} (p {v['p']:.3f})" for k, v in c.items() if k.endswith("r495")))
        out.append(f"   vs dpds: " + "; ".join(f"{k.split('~')[0]} {v['rho']:+.2f}" for k, v in c.items() if k.endswith("dpds495")))
    (HERE / "work/profile801.txt").write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
