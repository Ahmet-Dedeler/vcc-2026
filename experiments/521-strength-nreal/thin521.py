"""521 — is the cancer folds' weaker n_real prediction a POWER effect? The same S4 strength ranks against H1's own truth
rebuilt at fold power (480's gates: 150 / 100 / 250 cells, and 150 cells + 0.6 binomial thinning = the folds' depth).
    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/thin521.py
"""
import sys
sys.dont_write_bytecode = True
import json
import numpy as np
from pred521 import ROOT, WORK, evaluate, meanrank
G480 = ROOT / "experiments/480-nmae-instrument/work"
z = np.load(WORK / "pred521_h1.npz")
T = [str(t) for t in z["targets"]]
S4 = z["S4__S3_p_vip_thr4"] if "S4__S3_p_vip_thr4" in z.files else None
keys = [k for k in z.files if k.startswith("S4") or k.startswith("k562_t4")]
print("arrays:", keys)
feats = {"S4 thr4": z[[k for k in z.files if k.startswith("S4__S3") and k.endswith("thr4")][0]],
         "S4 thr-avg": z["S4_thr-avg"], "k562 thr4": z["k562_t4"]}
out = {}
for gname in ("h1_full_ka", "h1_c250", "h1_c150", "h1_c100", "h1_thin"):
    g = np.load(G480 / f"gate-{gname}.npz")
    Tg = [str(t) for t in g["targets"]]
    assert Tg == T, gname
    genes = g["genes"].astype(str)
    q = np.asarray(g["q"], float)
    n = np.array([int(((q[k] < 0.05) & (genes != t)).sum()) for k, t in enumerate(T)])
    print(f"{gname}: n_real median {np.median(n):.0f}, <20 {(n < 20).sum()}, 0 {(n == 0).sum()}")
    out[gname] = [evaluate(v, n, f"{nm} vs {gname}") for nm, v in feats.items()]
(WORK / "thin521.json").write_text(json.dumps(out, indent=1, default=float))
