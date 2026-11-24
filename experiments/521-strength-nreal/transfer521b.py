"""521 — transfer of the MODERATE shrinks (bottom q % by strength -> K 800 / 1,200 / 1,600; the smallest downside) with
the same percentile-utility curve as transfer521.py, plus the matched RANDOM selection of the same size (what the global
channel alone gives), so key-increment = selected - random is read per surface.
    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/transfer521b.py
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
n1 = P["nreal"]; p1 = pct(n1)
surf = {"H1|S4": (P["S4__S3_p_vip_thr4"], n1)}
for fold in ("k562", "hct116", "hek293t"):
    Tf, nf = fold_surface(fold)
    G = {b: strength(BANKS[b], Tf, 4.0) for b in ("k562", "hct116", "hek293t", "vip", "kolf2")}
    G["cd4rob"], G["cd4nde"] = cd4_feats(Tf)
    prod = [b for b in ("k562", "hct116", "hek293t", "vip") if b != fold and not (fold == "k562" and b == "vip")]
    surf[f"{fold}|all"] = (meanrank([G[b] for b in prod + ["kolf2", "cd4rob", "cd4nde"]]), nf)
out = {}
h = 0.08
for reg in ("1.0", "0.5", "0.35", "cmp"):
    R0 = T.rows[f"{reg}|base"]
    for K in (800, 1200, 1600):
        u = per_target_contrib(T.assemble(reg, "headcpm", "b", [K] * 85), R0)
        def curve(x, excl=None):
            w = np.exp(-0.5 * ((x[:, None] - p1[None, :]) / h) ** 2)
            if excl is not None:
                w[np.arange(x.size), excl] = 0.0
            return (w * u[None, :]).sum(1) / w.sum(1)
        line = []
        for s, (z, n) in surf.items():
            ok = np.isfinite(z); zz, nn = z[ok], n[ok]; pp = pct(nn)
            for q in (20, 50):
                k = int(round(q / 100 * zz.size))
                sel = np.argsort(zz, kind="stable")[:k]
                v_sel = float(curve(pp[sel], excl=sel if s.startswith("H1") else None).sum() * 85 / zz.size)
                v_rnd = float(curve(pp, excl=np.arange(zz.size) if s.startswith("H1") else None).mean() * k * 85 / zz.size)
                out[f"{reg}|K{K}|{s}|q{q}"] = dict(selected=v_sel, random=v_rnd, increment=v_sel - v_rnd)
                line.append(f"{s} q{q} {v_sel:+.4f} (rnd {v_rnd:+.4f})")
        print(f"{reg:4s} K {K:4d}: " + "  ".join(line))
(L.WORK / "transfer521b.json").write_text(json.dumps(out, indent=1))
