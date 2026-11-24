"""Does the incumbent emitter's delivered per-gene total deviate systematically from alloc(S, design profile)?"""
import sys, os
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
os.environ.setdefault("NUMBA_CACHE_DIR", str(HERE / "work" / "numba-cache"))
sys.path.insert(0, str(ROOT / "experiments/400-uncalled-knn")); sys.path.insert(0, str(HERE))
import run400 as H
R, E = H.R, H.E
b = H.builder389()
eff, genes, targets = H.native_cap15(b)
r, cc, emit_c = E.load_reference()
eff = R.map_axis(eff, genes, r["genes"])
S0 = np.load(HERE / "work/h1/stats-base395.npz")["sums"]
S1 = np.load(HERE / "work/h1/stats-bal_ident.npz")["sums"]
TS = 5e4
def lb(s): return np.log1p(TS * s / s.sum())
bins = [0, 0.5, 2, 5, 20, 100, 1e9]
acc = {k: [] for k in range(len(bins) - 1)}; accb = {k: [] for k in range(len(bins) - 1)}
for i in range(len(targets)):
    prof = b.profile_from_effect(eff[i], r["b_ctrl"], min_cpm=10.)
    tot = S0[i].sum(); exp = prof * tot
    d = lb(S0[i]) - np.log1p(TS * prof)          # delivered minus design, log-pseudobulk units
    db = lb(S1[i]) - np.log1p(TS * prof)
    for k in range(len(bins) - 1):
        m = (exp >= bins[k]) & (exp < bins[k + 1])
        acc[k].append(d[m]); accb[k].append(db[m])
print("expected-count bin    n_genes   mean(delivered-design) base395   bal_ident   (log-pseudobulk units)")
for k in range(len(bins) - 1):
    a = np.concatenate(acc[k]); c = np.concatenate(accb[k])
    print(f"[{bins[k]:>5},{bins[k+1]:>6})  {a.size:9d}   {a.mean():+.5f} (sd {a.std():.4f})   {c.mean():+.5f} (sd {c.std():.4f})")
