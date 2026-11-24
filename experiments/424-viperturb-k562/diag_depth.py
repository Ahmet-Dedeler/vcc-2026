"""424 diagnostic: is VIPerturb's weak agreement with K562-GWPS noise-limited? Per-target centred cosine (K562-GWPS gene
axis, target coordinates dropped) stratified by VIPerturb cells per target, for all shared targets and the panel."""
import sys
sys.dont_write_bytecode = True
import json
import numpy as np
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
def load(p):
    p = ROOT / p
    return (np.load(p / "delta.npy", mmap_mode="r"), [l for l in (p / "perts.txt").read_text().splitlines() if l],
            [l for l in (p / "genes.txt").read_text().splitlines() if l])
Dv, pv, gv = load("data/donors/viperturb_k562"); Dg, pg, gg = load("data/donors/k562_gwps")
cells = json.loads((ROOT / "data/donors/viperturb_k562/provenance.json").read_text())["cells_per_target"]
panel = [l.split(",")[0] for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
sh = sorted(set(pv) & set(pg)); ivp = {p: i for i, p in enumerate(pv)}; igp = {p: i for i, p in enumerate(pg)}
ivg = {g: i for i, g in enumerate(gv)}; igg = {g: i for i, g in enumerate(gg)}
genes = [g for g in gg if g in ivg and g not in set(panel)]   # drop the panel coordinates; own-target coordinate zeroed per row below
cv, cg = np.asarray([ivg[g] for g in genes]), np.asarray([igg[g] for g in genes])
V = np.stack([np.asarray(Dv[ivp[p]], dtype=np.float64)[cv] for p in sh]); G = np.stack([np.asarray(Dg[igp[p]], dtype=np.float64)[cg] for p in sh])
gpos = {g: i for i, g in enumerate(genes)}
for r, p in enumerate(sh):
    if p in gpos:
        V[r, gpos[p]] = 0.0; G[r, gpos[p]] = 0.0
V -= V.mean(0); G -= G.mean(0)
u = lambda M: M / np.linalg.norm(M, axis=1, keepdims=True)
cos = (u(V) * u(G)).sum(1)
n = np.asarray([cells[p] for p in sh]); pan = np.asarray([p in set(panel) for p in sh])
gn = np.linalg.norm(G, axis=1)
print(f"shared targets {len(sh)}, genes {len(genes)}; overall mean cos {cos.mean():+.4f}")
for lo, hi in ((20, 35), (35, 50), (50, 70), (70, 100), (100, 200)):
    m = (n >= lo) & (n < hi)
    print(f"  VIP cells [{lo:3d},{hi:3d}): n {m.sum():5d}  mean cos {cos[m].mean():+.4f}  median {np.median(cos[m]):+.4f}  (panel n {int((m & pan).sum())}: mean {cos[m & pan].mean() if (m & pan).any() else float('nan'):+.4f})")
q = np.quantile(gn, [0, .5, .8, .95, 1])
for lo, hi in zip(q[:-1], q[1:]):
    m = (gn >= lo) & (gn <= hi)
    print(f"  GWPS effect norm [{lo:5.2f},{hi:5.2f}]: n {m.sum():5d} mean cos {cos[m].mean():+.4f}  mean VIP cells {n[m].mean():.0f}")
top = np.argsort(-gn)[:300]
print(f"  top-300 strongest GWPS targets: mean cos {cos[top].mean():+.4f}; panel-254: mean cos {cos[pan].mean():+.4f}")
