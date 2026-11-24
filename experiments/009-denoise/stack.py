"""Do empirical-Bayes shrinkage and multi-donor averaging stack?

They attack different noise: EB removes per-gene sampling noise within a donor, multi-donor
averages independent noise across donors. If the noise sources are independent the gains
should compose.
"""
import sys, json
from pathlib import Path
import numpy as np, anndata as ad
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, TS

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

def prep(name, path, eb: bool, gene_norm: bool):
    """Return (full-axis effect matrix aligned to `genes`, index dict)."""
    d = Donor.load(path)
    m, s = d.decompose()
    dc = [j for j, g in enumerate(d.genes) if str(g) in gi]
    tc = [gi[str(d.genes[j])] for j in dc]
    dc, tc = np.asarray(dc), np.asarray(tc)
    if gene_norm:
        from vcc2026.donor_transfer import _norm_scale
        sc = _norm_scale(s); sn = s / sc
        s = (sn * (np.median(np.linalg.norm(s, 1 and 1, keepdims=False)) if False else 1.0)).astype(np.float32)
        s = sn.astype(np.float32)
    if eb:
        e = np.expm1(d.b_ctrl); pg = e / e.sum()
        meta = json.loads((Path(path) / "meta.json").read_text())
        n_cells = meta["cells_per_pert"]["median"]
        sigma2 = 1.0 / np.maximum(n_cells * 20000.0 * pg, 1e-3)
        V = s.var(0)
        tau2 = np.maximum(V - sigma2, 0.0)
        s = (s * (tau2 / (tau2 + sigma2 + 1e-12))).astype(np.float32)
    # rescale to a common effect magnitude so mixing weights mean the same thing
    s = s / (np.median(np.linalg.norm(s[:, dc], axis=1)) or 1.0)
    return s, dc, tc, {p: i for i, p in enumerate(d.perts)}

SPEC = {"k562_gwps": ("data/donors/k562_gwps", True),
        "orion_hct116": ("data/donors/orion_hct116", False)}

print("Transfer cosine on h1-test, K562 + Orion, with/without EB shrinkage\n")
results = {}
for eb in (False, True):
    D = {n: prep(n, p, eb, gn) for n, (p, gn) in SPEC.items()}
    sel = np.array([any(p in D[n][3] for n in D) for p in ps])
    for wk in (1.0, 2.0, 3.0, 1e9):
        w = {"k562_gwps": wk, "orion_hct116": 1.0}
        pred = np.zeros((int(sel.sum()), len(genes)))
        for r, p in enumerate(ps[sel]):
            have = [n for n in D if p in D[n][3]]
            tot = sum(w[n] for n in have) or 1.0
            for n in have:
                s, dc, tc, idx = D[n]
                pred[r, tc] += (w[n] / tot) * s[idx[p]][dc]
        c = cosm(pred, true[sel])
        label = "k562 only" if wk > 1e6 else f"{wk:g}:1"
        results[f"eb={eb} {label}"] = c
        print(f"  EB={str(eb):5}  weights {label:10}  cos {c:+.4f}")
best = max(results, key=results.get)
print(f"\nbest: {best} -> {results[best]:+.4f}")
base = results["eb=False k562 only"]
print(f"baseline (no EB, k562 only): {base:+.4f}   total gain {100*(results[best]-base)/base:+.1f}%")
Path("experiments/009-denoise/stack.json").write_text(json.dumps(results, indent=2))
