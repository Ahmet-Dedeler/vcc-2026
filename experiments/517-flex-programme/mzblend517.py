"""517: is a CLEAN4 blend a more ROBUST hook source than meanz_ne alone when the programme is attenuated? Hook alone on the
CLEAN4 block, gated; source = w z(meanz_ne, flipped at r) + (1 - w) z(CLEAN4) (CLEAN4 0 where uncovered), coverage matched to
Arc tau .05's; w 1 / 0.75 / 0.5; r 0 / 0.2 / 0.3 / 0.4; analytic, 3 seeds."""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import price517 as P
RC = P.RC
ESS = {l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt") if l.strip()}
zr = np.load(RC.REF, allow_pickle=False); genes_full = zr["genes"]; keep = np.asarray(zr["keep"], bool); T = zr["targets"].astype(str)
d = ROOT / "data/donors/h1_2025_full"; delta = np.load(d / "delta.npy"); perts = [l.strip() for l in open(d / "perts.txt")]
pos = {g: j for j, g in enumerate(l.strip() for l in open(d / "genes.txt"))}
col = np.array([pos[str(g)] for g in np.asarray(genes_full)[keep]])
rows = np.array([i for i, p in enumerate(perts) if p not in set(T)]); ne = np.array([perts[i] not in ESS for i in rows])
X = np.asarray(delta[rows][:, col], dtype=np.float64); sd = X.std(0); sd[sd == 0] = 1
z0 = (X[ne] / sd).mean(0); z0 = z0 / np.std(z0)
SRC, _ = P.load_sources(genes_full, keep)
zc4 = np.where(np.isfinite(SRC["c4"]), (0.5 - SRC["c4"]) / np.nanstd(SRC["c4"]), 0.0)
n_cov = int((np.isfinite(SRC["arc"]) & (np.abs(SRC["arc"] - 0.5) >= 0.05)).sum())
hc = np.where(np.isfinite(SRC["c4"]) & (np.abs(SRC["c4"] - 0.5) >= 0.05), np.where(SRC["c4"] > 0.5, -1.0, 1.0), 0.0)
c4 = RC.prior_signs_full(genes_full)[keep]; strong, _ = P.strong_set(T)
u = np.random.default_rng(51717).random(keep.sum())
res = {}
for seed in (0, 1, 2):
    ctx = P.Ctx("c4", seed, lambda i, t: c4); rng = np.random.default_rng(5183 + seed)
    for r in (0.0, 0.2, 0.3, 0.4):
        zr_ = np.where(u < r, -z0, z0)
        for w in (1.0, 0.75, 0.5):
            z = w * zr_ + (1 - w) * zc4
            a = np.abs(z); hs = np.where(a >= np.sort(a)[::-1][n_cov - 1], np.sign(z), 0.0)
            lfc, _ = ctx.hook(lambda i, t, hs=hs: hs if strong[i] else hc); sc, _ = ctx.score(lfc, rng, n_boot=200)
            res.setdefault((r, w), []).append((sc["nmae"], sc["reach"], sc["overall_equiv"]))
    del ctx
for (r, w), v in sorted(res.items()):
    a = np.array(v).mean(0)
    print(f"r {r:.1f} w {w:.2f}: nmae {a[0]:+.5f} reach {a[1]:+.5f} ov {a[2]:+.5f}")
