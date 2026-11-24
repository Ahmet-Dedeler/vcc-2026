"""450g: when the recipient's own line is ALREADY one of the donors (a K562 or HCT116 October context), what does the
same-line blend add over the donor mix? Truth = VIPerturb K562 (essential stratum, where the 48-cell truth has signal).
  mix3   = unit-mean of the three production donors' centred rows (K562-GWPS, HCT116, HEK293T): the same line at 1/3
  xline  = unit-mean of HCT116 + HEK293T only (what a line with no bank gets)
  same   = K562-GWPS alone
  blend(w) = unit(mix3) + w * unit(same)
Centred cosine to the truth (each matrix's mean over the evaluated targets removed), per target, mean."""
import runpy
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.argv = [sys.argv[0]]
src = open(HERE / "k562rep450.py").read().split('for sname, ts in strata.items():')[0]
g = {"__file__": str(HERE / "k562rep450.py"), "__name__": "k562rep450"}
exec(src, g)
B, strata, common_g, rng = g["B"], g["strata"], g["common_g"], np.random.default_rng(20260925)
ts = list(rng.choice(strata["essential"], size=min(1500, len(strata["essential"])), replace=False))
M = {}
for k, (delta, perts, genes) in B.items():
    pi = {p: i for i, p in enumerate(perts)}; gi = {x: i for i, x in enumerate(genes)}
    rows = np.array([pi[t] for t in ts]); cols = np.array([gi[x] for x in common_g]); order = np.argsort(rows)
    X = np.empty((len(ts), len(cols))); X[order] = np.asarray(delta[rows[order]], dtype=np.float64)[:, cols]
    own = {x: j for j, x in enumerate(common_g)}
    for r, t in enumerate(ts):
        if t in own:
            X[r, own[t]] = 0.0
    M[k] = X - X.mean(0)
u = lambda A: A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-12)
cos = lambda A, Bm: np.sum(u(A) * u(Bm), 1)
T = M["truth_viperturb"]
mix3 = u(u(M["k562_gwps"]) + u(M["hct116"]) + u(M["hek293t"]))
xline = u(u(M["hct116"]) + u(M["hek293t"]))
same = u(M["k562_gwps"])
print(f"essential targets {len(ts)}, genes {len(common_g)}; centred cosine to VIPerturb K562:")
for name, P in (("xline (no same-line bank)", xline), ("mix3 (same line at 1/3)", mix3), ("same (K562-GWPS alone)", same)):
    print(f"  {name:28s} {cos(P, T).mean():.4f}")
for w in (0.5, 1, 2, 4):
    print(f"  blend mix3 + {w:g}*same          {cos(u(mix3) + w * same, T).mean():.4f}")
