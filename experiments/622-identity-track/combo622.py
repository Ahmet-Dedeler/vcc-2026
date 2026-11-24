#!/usr/bin/env python3
"""622i — on the Broad-16 (12 working lines; reserved 4 never read): does the identity axis still add when a same-lab
3-line programme is available? Score = z(prog3 - .5) [+ lam z(I)] [+ mu z(C4 - .5)], gate-sign paired vs prog3 alone,
mean over 40 random triples per held-out line. Analysis only, < 4 GB."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import broad622 as B  # noqa: E402

RNG = np.random.default_rng(6229)


def main() -> None:
    gl, LFC, cl = B.load_matrix(B.P / "lfc_matrix.csv", 1)
    gf, FDR, cf = B.load_matrix(B.P / "fdr_matrix.csv", 1)
    gp, PB, cp = B.load_matrix(B.P / "pseudobulk_sum_by_ko.csv", 0)
    lines = sorted({ln for ln, _ in cl}); assert len(lines) == 12 and not (set(lines) & B.RESERVED)
    c4 = np.load(B.ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    nm = list(c4["bank_names"].astype(str))
    Pc = c4["P_banks"].astype(float)[[nm.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    cmap = dict(zip(c4["genes"].astype(str), Pc)); C4 = np.array([cmap.get(g, np.nan) for g in gl])
    zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
    all_t = set(k for _, k in cl)
    gpos = {g: i for i, g in enumerate(gp)}
    cpm = {}
    for ln in lines:
        cc = [j for j, (l2, k) in enumerate(cp) if l2 == ln and k not in all_t]
        s = np.nansum(PB[:, cc], 1); v = 1e6 * s / np.nansum(s); v[np.all(~np.isfinite(PB[:, cc]), 1)] = np.nan; cpm[ln] = v
    rowmap = np.array([gpos.get(g, -1) for g in gl])
    gated = (FDR < 0.05) & np.isfinite(LFC) & (LFC != 0); dn = (LFC < 0) & gated
    colline = np.array([ln for ln, _ in cl])
    pdn = {ln: dn[:, colline == ln].sum(1) for ln in lines}; pg = {ln: gated[:, colline == ln].sum(1) for ln in lines}
    out = {}
    for L in lines:
        others = [o for o in lines if o != L]
        Ls = np.log1p(cpm[L]); Lr = np.nanmean([np.log1p(cpm[o]) for o in others], 0)
        I = np.where(rowmap >= 0, (Ls - Lr)[np.maximum(rowmap, 0)], np.nan)
        I = np.where(np.where(rowmap >= 0, cpm[L][np.maximum(rowmap, 0)], np.nan) > 5, I, np.nan)
        zi = (I - np.nanmean(I)) / np.nanstd(I)
        jl = np.flatnonzero(colline == L); tl = [cl[j][1] for j in jl]
        gate = gated[:, jl] & (gl[:, None] != np.array(tl)[None, :]); tdn = LFC[:, jl] < 0
        acc = {k: [] for k in ("prog3", "prog3+0.5I", "prog3+1I", "prog3+0.5C4", "prog3+0.5C4+0.5I")}
        for _ in range(40):
            tri = list(RNG.choice(others, 3, replace=False))
            sd = sum(pdn[o] for o in tri); sg = sum(pg[o] for o in tri)
            pr = (sd + 0.5) / (sg + 1.0)
            zp = np.where(sg >= 3, (pr - 0.5) / np.nanstd(pr - 0.5), np.nan)
            ok = np.isfinite(zp) & np.isfinite(zi) & np.isfinite(zc)
            m = gate & ok[:, None]
            base = None
            for k, sc in (("prog3", zp), ("prog3+0.5I", zp + 0.5 * zi), ("prog3+1I", zp + 1.0 * zi),
                          ("prog3+0.5C4", zp + 0.5 * zc), ("prog3+0.5C4+0.5I", zp + 0.5 * zc + 0.5 * zi)):
                a = float((((sc > 0)[:, None] == tdn) & m).sum() / max(m.sum(), 1))
                acc[k].append(a)
        res = {k: float(np.mean(v)) for k, v in acc.items()}
        res.update({f"d_{k}": res[k] - res["prog3"] for k in acc if k != "prog3"})
        out[L] = res
        print(f"{L:9s} prog3 {res['prog3']:.3f} | +0.5I {res['d_prog3+0.5I']:+.4f}  +1I {res['d_prog3+1I']:+.4f}  +0.5C4 {res['d_prog3+0.5C4']:+.4f}  +0.5C4+0.5I {res['d_prog3+0.5C4+0.5I']:+.4f}", flush=True)
    summ = {k: float(np.mean([out[L][k] for L in lines])) for k in out[lines[0]] if k.startswith("d_")}
    pos = {k: int(sum(out[L][k] > 0 for L in lines)) for k in summ}
    out["_summary"] = {"mean": summ, "lines_positive": pos}
    print("\nmean Δ vs prog3 alone:", {k: round(v, 4) for k, v in summ.items()}, "\nlines > 0:", pos)
    (HERE / "work/combo622.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
