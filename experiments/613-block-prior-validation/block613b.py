#!/usr/bin/env python3
"""613b — DESCRIPTIVE refinement of block613 (not used by any pre-registered rule): the builder's `clean4q80` applies a
block prior only to block genes with |E0| <= that target's block 80th percentile (scripts/make_submission_pm_cd4.py:622),
and the strongest keep their delivered sign. This re-reads identity (v3 side, λ 0.5) and 458's programme vs CLEAN4 on
those APPLIED block cells only, then converts each target's change in correct-sign calls into a fid estimate:
Δfid ≈ Δcorrect / 2,600 calls (n_pred >= n_real), averaged over targets.
E0 proxy: H1 = lib521's design key (the production key); folds = the 3-donor consensus key (cons0, as ds612/decomp639).
Approximations, stated: v3's confidence map and the builder's shrink/τ can leave some applied cells unflipped (458d's
dose was 0.43). This is an upper bound on the exposure at full dose.
Output: work/block613b.out.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import block613 as B  # noqa: E402  (puts 639 / 636 / 521 / 612 on the path)
import decomp639 as DC  # noqa: E402
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402
import lib521 as L  # noqa: E402

NCALL = 2600


def e0_for(I, s):
    D = I.S[s]
    T, G = D["T"], D["G"]
    if s == "h1":
        st = L.load_state(0)
        pos = {g: i for i, g in enumerate(st["genes"])}
        cols = np.array([pos.get(g, -1) for g in G])
        E = np.zeros((len(T), len(G)))
        E[:, cols >= 0] = np.abs(st["key"][:, cols[cols >= 0]])
        return E
    import ds612 as DS
    Sx = DS.F.load_surface(DC.FOLD[s])
    genes_f, targets_f = [str(g) for g in Sx["genes"]], list(Sx["targets"])
    cons0 = np.nan_to_num(np.nanmean(np.stack([DS.load_donor_rows(d, targets_f, np.array(genes_f), pc_cpm=1.0)[0]
                                               for d in Sx["donors"]]), 0))
    tpos, gpos = {t: i for i, t in enumerate(targets_f)}, {g: i for i, g in enumerate(genes_f)}
    ti, gi = np.array([tpos.get(t, -1) for t in T]), np.array([gpos.get(g, -1) for g in G])
    E = np.zeros((len(T), len(G)))
    for a, t in enumerate(ti):
        if t >= 0:
            E[a, gi >= 0] = np.abs(cons0[t, gi[gi >= 0]])
    return E


def main():
    I = Inst()
    prof = DC.flex_profiles()
    z458 = np.load(B.A458, allow_pickle=True)
    g458, p458 = z458["genes"].astype(str), np.asarray(z458["pdown_A"], float)
    lines = []
    for s in ("h1",) + B.CANCER:
        D = I.S[s]
        T, G, lfc, gate, tested = D["T"], D["G"], D["lfc"], D["gate"], D["tested"]
        q = np.load(G480 / SURF[s], allow_pickle=True)["q"].astype(float)
        null = tested & np.isfinite(q) & (q >= 0.05) & (lfc != 0)
        if s == "h1":
            own, refs = prof["h1"], [prof[k] for k in ("ctxA", "ctxB", "ctxC")]
        else:
            zz = np.load(G480 / SURF[s], allow_pickle=True)
            own = dict(zip(zz["genes"].astype(str), zz["ctrl_cpm"].astype(float)))
            refs = [prof[k] for k in ("ctxA", "ctxB", "ctxC", "h1")]
        Iv = DC.identity(own, refs, list(G))
        zi = (Iv - np.nanmean(Iv)) / np.nanstd(Iv)
        C4, _ = clean4_lolo(I, s)
        zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
        sc = zc + B.LAM * zi
        _, block = DC.callset_for(I, s)
        E = e0_for(I, s)
        apply = np.zeros_like(block)
        for i in range(len(T)):
            b = block[i]
            if b.any():
                apply[i] = b & (E[i] <= np.quantile(E[i, b], 0.8))
        dn = lfc < 0
        P458 = I.gene_vec(s, g458, p458)
        lines.append(f"== {s}: applied block cells per target (median) {int(np.median(apply.sum(1)))} of {int(np.median(block.sum(1)))}")
        for lab, src, dfn in (("identity", sc > 0, np.isfinite(zc) & np.isfinite(zi)), ("458", P458 > 0.5, np.isfinite(zc) & np.isfinite(P458))):
            dc_t = np.zeros(len(T))
            parts = []
            for rlab, M in (("gate", gate & apply), ("null", null & apply)):
                m = M & dfn[None, :]
                sgood = ((src[None, :] == dn) & m).sum(1)
                cgood = (((zc > 0)[None, :] == dn) & m).sum(1)
                dc_t += sgood - cgood
                n = m.sum()
                parts.append(f"{rlab} {float((sgood.sum() - cgood.sum()) / max(n, 1)):+.4f} (n/target {n / len(T):.0f})")
            lines.append(f"  {lab:8s} " + "  ".join(parts) + f"  -> Δfid ≈ {float(np.mean(dc_t)) / NCALL:+.4f} raw "
                         f"({float(np.mean(dc_t)) / NCALL * 3.32482 / 6:+.4f} overall)")
    print("\n".join(lines))
    (HERE / "work/block613b.out").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
