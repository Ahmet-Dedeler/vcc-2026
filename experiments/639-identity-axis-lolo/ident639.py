#!/usr/bin/env python3
"""639 — does 441's "identity genes go down" signal transfer leave-one-line-out? (council 600, DE wall, portfolio item 2)

441 PS6 (h1-test, same chemistry): the identity axis I = log ctrl_line − mean log ctrl_others correlates −0.236 with H1's
mean response; combined with CLEAN4 (z(CLEAN4) + λ z(−I), λ 1–2) it lifted gate-block sign accuracy 0.580 → 0.592. For
D/E/F it is computable on the day, within Flex, from the controls alone (D/E/F vs A/B/C and each other). Never built,
because H1 alone priced it at ~+0.001 on the board.

Here: on 636's four gate surfaces, each with its own control CPM profile (480's gate files: `ctrl_cpm`), plus the A/B/C
challenge controls (Flex), the identity of surface S = log(ctrl_S + 1) − mean over the OTHER profiles of log(ctrl + 1),
on genes all profiles measure. Score = z(P_C4 − 0.5) + λ z(I) (positive → predict down), with CLEAN4-LOLO; identity alone
fills CLEAN4's coverage gaps. Reported: paired Δ vs CLEAN4 on CLEAN4's genes (λ sweep) and the full-coverage reading
(identity-filled vs all-down-filled). Analysis only, < 2 GB.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402

WORK = HERE / "work"
WORK.mkdir(exist_ok=True)


def abc_profiles():
    import anndata as ad
    genes = __import__("pandas").read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
    out = {}
    for c in "ABC":
        a = ad.read_h5ad(ROOT / f"data/challenge/context_{c}.h5ad", backed="r")
        X = a.X[:]
        lib = np.asarray(X.sum(1)).ravel()
        import scipy.sparse as sp
        cpm = np.asarray((sp.diags(1e6 / lib) @ X).mean(0)).ravel()
        out[f"ctx{c}"] = dict(zip(genes, cpm))
        del X
    return out


def main():
    I = Inst()
    prof = {}
    for s, f in SURF.items():
        z = np.load(G480 / f, allow_pickle=True)
        prof[s] = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(float)))
    prof.update(abc_profiles())
    common = set.intersection(*(set(p) for p in prof.values()))
    print(f"profiles {list(prof)}; genes measured by all: {len(common)}", flush=True)
    cg = sorted(common)
    L = {k: np.log1p(np.array([prof[k][g] for g in cg])) for k in prof}
    rows = []
    for s in I.S:
        others = [k for k in prof if k != s]
        ident = L[s] - np.mean([L[k] for k in others], 0)
        Iv = I.gene_vec(s, cg, ident)                    # NaN where the surface's gene is not in the common set
        C4, _ = clean4_lolo(I, s)
        zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
        zi = (Iv - np.nanmean(Iv)) / np.nanstd(Iv)
        D = I.S[s]
        gate, dn = D["gate"], D["lfc"] < 0
        def acc(score, mask):
            m = gate & mask[None, :] & np.isfinite(score)[None, :] & (score != 0)[None, :]
            ok = ((score > 0)[None, :] == dn) & m
            return float(ok.sum() / max(m.sum(), 1)), int(m.sum())
        both = np.isfinite(zc) & np.isfinite(zi)
        r = dict(surface=s, n_common=int(np.isfinite(Iv).sum()))
        r["corr_identity_vs_C4"] = float(np.corrcoef(zc[both], zi[both])[0, 1])
        r["C4_on_both"] = acc(zc, both)
        r["identity_alone_on_both"] = acc(zi, both)
        for lam in (0.25, 0.5, 1.0, 2.0):
            r[f"C4+{lam}I_on_both"] = acc(zc + lam * zi, both)
        gaps = ~np.isfinite(zc) & np.isfinite(zi)
        r["identity_on_C4_gaps"] = acc(zi, gaps)
        r["alldown_on_C4_gaps"] = acc(np.ones_like(zi), gaps)
        rows.append(r)
        print(json.dumps(r), flush=True)
    (WORK / "ident639.json").write_text(json.dumps(rows, indent=1))
    print(f"\n{'surface':8s} {'C4':>13s} {'I alone':>13s} {'C4+.5I':>13s} {'C4+1I':>13s} {'C4+2I':>13s} | gaps: {'I':>13s} {'all-down':>13s}  r(C4,I)")
    for r in rows:
        f = lambda k: f"{r[k][0]:.4f} ({r[k][1]})"
        print(f"{r['surface']:8s} {f('C4_on_both'):>13s} {f('identity_alone_on_both'):>13s} {f('C4+0.5I_on_both'):>13s} "
              f"{f('C4+1.0I_on_both'):>13s} {f('C4+2.0I_on_both'):>13s} | {f('identity_on_C4_gaps'):>13s} "
              f"{f('alldown_on_C4_gaps'):>13s}  {r['corr_identity_vs_C4']:+.3f}")


if __name__ == "__main__":
    sys.exit(main())
