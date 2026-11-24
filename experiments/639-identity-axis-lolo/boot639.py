#!/usr/bin/env python3
"""639b - target-bootstrap CIs for CLEAN4 + lambda*identity vs CLEAN4 (paired, genes where both are defined), lambda FIXED
at 1.0 (441 PS6's natural equal-weight choice; 0.5 and 2 reported as sensitivity), and the full-coverage reading:
[combined where CLEAN4 exists, identity alone in CLEAN4's gaps] vs [CLEAN4 where it exists, nothing elsewhere]."""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import G480, SURF, Inst, clean4_lolo
import ident639 as M
I = Inst()
prof = {}
for s, f in SURF.items():
    z = np.load(G480 / f, allow_pickle=True); prof[s] = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(float)))
prof.update(M.abc_profiles())
cg = sorted(set.intersection(*(set(p) for p in prof.values())))
L = {k: np.log1p(np.array([prof[k][g] for g in cg])) for k in prof}
rng = np.random.default_rng(639); NB = 4000
rows = []
for s in I.S:
    ident = L[s] - np.mean([L[k] for k in prof if k != s], 0)
    Iv = I.gene_vec(s, cg, ident); C4, _ = clean4_lolo(I, s)
    zc = (C4 - 0.5) / np.nanstd(C4 - 0.5); zi = (Iv - np.nanmean(Iv)) / np.nanstd(Iv)
    D = I.S[s]; gate, dn = D["gate"], D["lfc"] < 0
    both = np.isfinite(zc) & np.isfinite(zi)
    m = gate & both[None, :]
    base_ok = ((zc > 0)[None, :] == dn) & m
    for lam in (0.5, 1.0, 2.0):
        sc = zc + lam * zi
        ok = ((sc > 0)[None, :] == dn) & m
        a, b, n = ok.sum(1).astype(float), base_ok.sum(1).astype(float), m.sum(1).astype(float)
        k = n > 0; a, b, n = a[k], b[k], n[k]
        idx = rng.integers(0, len(n), size=(NB, len(n)))
        d = (a[idx].sum(1) - b[idx].sum(1)) / n[idx].sum(1)
        rows.append(dict(surface=s, lam=lam, reading="paired on CLEAN4+identity genes", d=float((a.sum()-b.sum())/n.sum()),
                         ci=[float(np.quantile(d, .025)), float(np.quantile(d, .975))], p_le0=float((d <= 0).mean()),
                         n_targets=int(k.sum()), n_pairs=int(n.sum())))
    # full-coverage: combined (lam 1) on both + identity alone in CLEAN4 gaps, vs CLEAN4 alone where defined (no call elsewhere)
    sc = np.where(np.isfinite(zc), zc + np.nan_to_num(zi), zi)
    covA = gate & np.isfinite(sc)[None, :]
    covB = gate & np.isfinite(zc)[None, :]
    okA = ((sc > 0)[None, :] == dn) & covA
    okB = ((zc > 0)[None, :] == dn) & covB
    # net correct-minus-wrong calls per gate gene (what an MDE-sized hook earns, sign-wise), per target
    netA = (2 * okA.sum(1) - covA.sum(1)).astype(float); netB = (2 * okB.sum(1) - covB.sum(1)).astype(float)
    ng = gate.sum(1).astype(float); k = ng > 0
    idx = rng.integers(0, int(k.sum()), size=(NB, int(k.sum())))
    dd = (netA[k][idx].sum(1) - netB[k][idx].sum(1)) / ng[k][idx].sum(1)
    rows.append(dict(surface=s, lam=1.0, reading="full coverage: net (right-wrong)/gate genes, combined+gaps vs CLEAN4 only",
                     d=float((netA[k].sum() - netB[k].sum()) / ng[k].sum()), ci=[float(np.quantile(dd, .025)), float(np.quantile(dd, .975))],
                     p_le0=float((dd <= 0).mean()), cov_combined=float(covA.sum() / gate.sum()), cov_c4=float(covB.sum() / gate.sum())))
for r in rows:
    extra = f"cov {r['cov_c4']:.2f}->{r['cov_combined']:.2f}" if 'cov_combined' in r else f"targets {r['n_targets']} pairs {r['n_pairs']}"
    print(f"{r['surface']:8s} lam {r['lam']:.1f} {r['reading'][:44]:44s} d {r['d']:+.4f} [{r['ci'][0]:+.4f}, {r['ci'][1]:+.4f}] P(<=0) {r['p_le0']:.3f}  {extra}")
(HERE / "work/boot639.json").write_text(json.dumps(rows, indent=1))
