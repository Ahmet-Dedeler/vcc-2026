"""450b: the same comparison on the genes BOTH priors cover, and per target (bootstrap CI). Imports sameline450's setup."""
import runpy
import sys
from pathlib import Path

import numpy as np

sys.argv = [sys.argv[0]]
g = runpy.run_path(str(Path(__file__).resolve().parent / "sameline450.py"), run_name="not_main")
P, lfc, gate, block, adj, truth, T = (g[k] for k in ("P", "lfc", "gate", "block", "adj", "truth", "T"))
both = np.isfinite(P["clean4"]) & np.isfinite(P["pan_pdown"])
L = np.where(np.isfinite(lfc), lfc, 0.0)
S = L.sum(0)
loo = np.sign(S[None, :] - L)
for name, sub in (("gate_block", gate & block), ("gate_all", gate), ("block_all_adj", adj & block)):
    m = sub & both[None, :]
    r = {k: float(((np.broadcast_to(P[k], lfc.shape) == truth) & m).sum() / m.sum()) for k in ("clean4", "pan_pdown")}
    r["oracle"] = float(((loo == truth) & m).sum() / m.sum())
    print(f"{name:14s} common pairs {int(m.sum()):6d}  " + "  ".join(f"{k} {v:.4f}" for k, v in r.items()))
acc = {k: [] for k in ("clean4", "pan_pdown", "oracle")}
for i in range(len(T)):
    m = gate[i] & block & both
    if m.sum() < 10:
        continue
    acc["clean4"].append(np.mean(P["clean4"][m] == truth[i][m]))
    acc["pan_pdown"].append(np.mean(P["pan_pdown"][m] == truth[i][m]))
    acc["oracle"].append(np.mean(loo[i][m] == truth[i][m]))
a = {k: np.array(v) for k, v in acc.items()}
d = a["pan_pdown"] - a["clean4"]
rng = np.random.default_rng(0)
bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(4000)]
print(f"per target (gate-block >= 10 common genes): n {len(d)}; means " + ", ".join(f"{k} {v.mean():.4f}" for k, v in a.items())
      + f"; pan - clean4 {d.mean():+.4f} [95% {np.percentile(bs, 2.5):+.4f}, {np.percentile(bs, 97.5):+.4f}]; pan better on "
      f"{int((d > 0).sum())} / {len(d)}")
