"""Target-sampling bootstrap (paired over the scored targets, seed-averaged per-target deltas) for nmae."""
import sys, numpy as np, polars as pl
from pathlib import Path
HERE = Path(__file__).resolve().parent
W = HERE / "work/h1"
sfx = sys.argv[2] if len(sys.argv) > 2 else ""
rng = np.random.default_rng(423)
def per(tag):
    d = pl.read_csv(W / f"per-target-{tag}.csv").filter(pl.col("metric") == "de_wilcoxon_lfc_nmae")
    return dict(zip(d["perturbation"].to_list(), d["value"].to_list()))
for v in sys.argv[1].split(","):
    deltas = []
    for s in (0, 1, 2):
        f = W / f"per-target-{v}_s{s}{sfx}.csv"
        if not f.exists():
            continue
        b, x = per(f"base_s{s}{sfx}"), per(f"{v}_s{s}{sfx}")
        ks = sorted(set(b) & set(x))
        deltas.append(np.array([x[k] - b[k] for k in ks]))
    D = np.mean(deltas, axis=0)
    bs = np.array([rng.choice(D, D.size).mean() for _ in range(20000)])
    print(f"{v}{sfx}: seeds {len(deltas)}, targets {D.size}, mean {D.mean():+.5f}, 95% CI [{np.percentile(bs, 2.5):+.5f}, {np.percentile(bs, 97.5):+.5f}], "
          f"p(>=0) {np.mean(bs >= 0):.4f}; targets better {int((D < -1e-6).sum())}, worse {int((D > 1e-6).sum())}, "
          f"top-3 targets carry {np.sort(D)[:3].sum() / D.sum():.2f} of the gain")
