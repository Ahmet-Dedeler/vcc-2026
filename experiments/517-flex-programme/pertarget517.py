"""517: per-target breadth of the hook's nmae gain on H1 (analytic, seed 0): meanz_ne vs P_down vs CLEAN4 sources, gated,
CLEAN4 block. How many of the nmae-scored targets improve, and how concentrated is the gain (top-5 share)?"""
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
z = (X[ne] / sd).mean(0)
SRC, _ = P.load_sources(genes_full, keep)
n_cov = int((np.isfinite(SRC["arc"]) & (np.abs(SRC["arc"] - 0.5) >= 0.05)).sum())
a = np.abs(z); hmz = np.where(a >= np.sort(a)[::-1][n_cov - 1], np.sign(z), 0.0)
harc = np.where(np.isfinite(SRC["arc"]) & (np.abs(SRC["arc"] - 0.5) >= 0.05), np.where(SRC["arc"] > 0.5, -1.0, 1.0), 0.0)
hc = np.where(np.isfinite(SRC["c4"]) & (np.abs(SRC["c4"] - 0.5) >= 0.05), np.where(SRC["c4"] > 0.5, -1.0, 1.0), 0.0)
c4 = RC.prior_signs_full(genes_full)[keep]; strong, _ = P.strong_set(T)
ctx = P.Ctx("c4", 0, lambda i, t: c4); rng = np.random.default_rng(0)
for name, hs in (("meanz_ne", hmz), ("P_down", harc), ("CLEAN4", hc)):
    lfc, _ = ctx.hook(lambda i, t, hs=hs: hs if strong[i] else hc)
    _, D = ctx.score(lfc, rng, n_boot=10)
    dn = D["nmae"]; ok = np.isfinite(dn)
    v = dn[ok]; order = np.sort(v)
    ess = np.array([t in ESS for t in T])[ok]
    print(f"{name:9s} nmae-scored targets {ok.sum()}: improved {int((v < 0).sum())}, worse {int((v > 0).sum())}; mean {v.mean():+.5f}, "
          f"median {np.median(v):+.5f}; top-5 share of the gain {order[:5].sum() / v.sum():.2f}; non-essential mean {v[~ess].mean():+.5f} "
          f"(n {int((~ess).sum())}), essential {v[ess].mean():+.5f} (n {int(ess.sum())})")
