import numpy as np, polars as pl, sys
from pathlib import Path
W = Path(__file__).resolve().parent / "work/h1"
ROOT = Path(__file__).resolve().parents[2]
ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
keep = ref["keep"]; rq = ref["q"]; rl = ref["lfc"]
A = np.load(W / "stats-base395.npz"); B = np.load(W / f"stats-{sys.argv[1]}.npz")
pa = pl.read_csv(W / "per-target-base395.csv"); pb = pl.read_csv(W / f"per-target-{sys.argv[1]}.csv")
ra = pa.filter(pl.col("metric") == "de_wilcoxon_direction_reach_raw").sort("perturbation")
rb = pb.filter(pl.col("metric") == "de_wilcoxon_direction_reach_raw").sort("perturbation")
d = rb["value"].to_numpy() - ra["value"].to_numpy()
print("per-target reach delta: mean %+.4f  n<0 %d  n>0 %d  n=0 %d  q05 %+.3f q95 %+.3f" % (np.nanmean(d), (d < 0).sum(), (d > 0).sum(), (d == 0).sum(), np.nanquantile(d, .05), np.nanquantile(d, .95)))
# sign agreement of delivered lfc between the two arms, on the reference's significant genes, and p-value shifts
qa, qb, la, lb = A["q"], B["q"], A["lfc"], B["lfc"]
flip = []; sigchg = []
for i in range(qa.shape[0]):
    m = (rq[i] < 0.05) & np.isfinite(rl[i])
    flip.append(np.mean(np.sign(la[i][m]) != np.sign(lb[i][m])) if m.any() else np.nan)
    sigchg.append(np.mean((qa[i] < 0.05) != (qb[i] < 0.05)))
print("frac of reference-DE genes whose delivered LFC sign differs: mean %.4f max %.4f" % (np.nanmean(flip), np.nanmax(flip)))
print("frac of genes whose call (q<0.05) differs: mean %.4f max %.4f" % (np.mean(sigchg), np.max(sigchg)))
