"""521 — transfer of 'bottom q % by strength -> K' from H1 to a surface with a different tail: each selected target's
utility is read off H1's utility-vs-truth-percentile curve (LOTO Gaussian kernel regression on H1's per-target utilities,
bandwidth 0.08 in percentile units), at that surface's selected targets' truth percentiles. H1 itself is scored LOTO
(its own target excluded from the curve), so H1's row is an out-of-sample check of the curve against its realized value.
    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/transfer521.py
"""
import sys
sys.dont_write_bytecode = True
import json
import numpy as np
from scipy.stats import rankdata
import lib521 as L
from policy521 import Table, per_target_contrib, rank_policy, value
from pred521 import BANKS, cd4_feats, fold_surface, meanrank, strength

def pct(n):
    return (rankdata(n.astype(float), method="average") - 1) / (n.size - 1)

T = Table(0)
P = np.load(L.WORK / "pred521_h1.npz")
n1 = P["nreal"]
p1 = pct(n1)
surf = {"H1|S4": (P["S4__S3_p_vip_thr4"], n1), "H1|S4avg": (P["S4_thr-avg"], n1)}
for fold in ("k562", "hct116", "hek293t"):
    Tf, nf = fold_surface(fold)
    G = {b: strength(BANKS[b], Tf, 4.0) for b in ("k562", "hct116", "hek293t", "vip", "kolf2")}
    G["cd4rob"], G["cd4nde"] = cd4_feats(Tf)
    prod = [b for b in ("k562", "hct116", "hek293t", "vip") if b != fold and not (fold == "k562" and b == "vip")]
    surf[f"{fold}|prod"] = (meanrank([G[b] for b in prod]), nf)
    surf[f"{fold}|all"] = (meanrank([G[b] for b in prod + ["kolf2", "cd4rob", "cd4nde"]]), nf)
out = {}
h = 0.08
for reg in ("1.0", "cmp"):
    for var in ("b", "a"):
        R0 = T.rows[f"{reg}|base"]
        for K in (0, 10, 25, 100, 400):
            u = per_target_contrib(T.assemble(reg, "headcpm", var, [K] * 85), R0)
            def curve(x, excl=None):
                w = np.exp(-0.5 * ((x[:, None] - p1[None, :]) / h) ** 2)
                if excl is not None:
                    w[np.arange(x.size), excl] = 0.0
                return (w * u[None, :]).sum(1) / w.sum(1)
            row = {}
            for s, (z, n) in surf.items():
                ok = np.isfinite(z)
                zz, nn = z[ok], n[ok]
                pp = pct(nn)
                for q in (5, 10, 20):
                    k = max(1, int(round(q / 100 * zz.size)))
                    sel = np.argsort(zz, kind="stable")[:k]
                    uu = curve(pp[sel], excl=sel if s.startswith("H1") else None)
                    val = float(uu.sum() * 85 / zz.size)      # rescale the per-target sum to a panel of H1's size
                    row[f"{s}|q{q}"] = val
                    if s.startswith("H1"):
                        realized = value(T.assemble(reg, "headcpm", var, rank_policy(z, q, K)), R0)
                        row[f"{s}|q{q}|realized"] = realized
            out[f"{reg}|{var}|K{K}"] = row
            print(f"{reg:4s} {var} K {K:4d}: " + "  ".join(
                f"{s} q5 {row[f'{s}|q5']:+.4f} q10 {row[f'{s}|q10']:+.4f}" + (f" (H1 realized q10 {row[f'{s}|q10|realized']:+.4f})" if s.startswith('H1') else '')
                for s in surf))
(L.WORK / "transfer521.json").write_text(json.dumps(out, indent=1))
