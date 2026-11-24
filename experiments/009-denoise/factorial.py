"""Clean factorial: {raw, gene_norm, EB shrinkage} x {K562 only, K562+Orion 2:1}.

The first stack attempt applied gene_norm AND EB together and saw no EB gain. Both
downweight noisy genes, so they do overlapping work -- this separates them.
"""
import sys, json
from pathlib import Path
import numpy as np, anndata as ad
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, TS, _norm_scale

a = ad.read_h5ad("data/harness/h1-test.h5ad"); X = a.X.tocsr()
lab = a.obs["target_gene"].astype(str).to_numpy(); genes = np.asarray(a.var_names)
Pc = np.asarray(X[lab == "non-targeting"].sum(0)).ravel().astype(np.float64)
b = np.log1p(TS * Pc / Pc.sum())
ps, ds = [], []
for g in sorted(set(lab) - {"non-targeting"}):
    P = np.asarray(X[lab == g].sum(0)).ravel().astype(np.float64)
    ds.append(np.log1p(TS * P / P.sum()) - b); ps.append(g)
true = np.vstack(ds); ps = np.asarray(ps)
gi = {g: i for i, g in enumerate(genes)}

def cosm(P, Q):
    Pc_ = P - P.mean(0, keepdims=True); Qc = Q - Q.mean(0, keepdims=True)
    na, nb = np.linalg.norm(Pc_, axis=1), np.linalg.norm(Qc, axis=1)
    ok = (na > 0) & (nb > 0)
    return float(np.mean(np.sum(Pc_[ok]*Qc[ok], axis=1)/(na[ok]*nb[ok])))

def prep(path, treat):
    d = Donor.load(path)
    m, s = d.decompose()
    dc = np.asarray([j for j, g in enumerate(d.genes) if str(g) in gi])
    tc = np.asarray([gi[str(d.genes[j])] for j in dc])
    if treat == "gene_norm":
        s = (s / _norm_scale(s)).astype(np.float32)
    elif treat == "eb":
        e = np.expm1(d.b_ctrl); pg = e / e.sum()
        n_cells = json.loads((Path(path) / "meta.json").read_text())["cells_per_pert"]["median"]
        sigma2 = 1.0 / np.maximum(n_cells * 20000.0 * pg, 1e-3)
        V = s.var(0); tau2 = np.maximum(V - sigma2, 0.0)
        s = (s * (tau2 / (tau2 + sigma2 + 1e-12))).astype(np.float32)
    s = s / (np.median(np.linalg.norm(s[:, dc], axis=1)) or 1.0)
    return s, dc, tc, {p: i for i, p in enumerate(d.perts)}

PATHS = {"k562_gwps": "data/donors/k562_gwps", "orion_hct116": "data/donors/orion_hct116"}
print(f"{'treatment':14} {'K562 only':>11} {'K562+Orion 2:1':>16}")
out = {}
for treat in ("raw", "gene_norm", "eb"):
    D = {n: prep(p, treat) for n, p in PATHS.items()}
    sel = np.array([any(p in D[n][3] for n in D) for p in ps])
    row = []
    for w in ({"k562_gwps": 1.0, "orion_hct116": 0.0}, {"k562_gwps": 2.0, "orion_hct116": 1.0}):
        pred = np.zeros((int(sel.sum()), len(genes)))
        for r, p in enumerate(ps[sel]):
            have = [n for n in D if p in D[n][3] and w[n] > 0]
            tot = sum(w[n] for n in have) or 1.0
            for n in have:
                s, dc, tc, idx = D[n]
                pred[r, tc] += (w[n] / tot) * s[idx[p]][dc]
        row.append(cosm(pred, true[sel]))
    out[treat] = row
    print(f"{treat:14} {row[0]:11.4f} {row[1]:16.4f}")
base = out["raw"][0]
best = max((v, t, i) for t, r in out.items() for i, v in enumerate(r))
lbl = ["K562 only", "K562+Orion 2:1"][best[2]]
print(f"\nbaseline (raw, K562 only) {base:.4f}")
print(f"best: {best[1]} + {lbl} -> {best[0]:.4f}  ({100*(best[0]-base)/base:+.1f}%)")
Path("experiments/009-denoise/factorial.json").write_text(json.dumps(out, indent=2))
