#!/usr/bin/env python3
"""622 — the identity sign source (639) in its BOARD / DAY-OF form: reference set, stability, halves.

639: score = z(P_C4 - .5) + lam * z(I), I = log1p(ctrl CPM of line S) - mean over OTHER profiles of log1p(ctrl CPM),
positive -> predict down; beats CLEAN4 LOLO on 3/4 surfaces with 6 mixed-chemistry reference profiles.
On the board the references can only be the other Arc Flex contexts (A/B/C, optionally H1 Flex); on Oct 22, D/E/F vs
A/B/C + the other two. This script asks:
  (a) h1 with FLEX-ONLY references {A, B, C} (exactly the day-of construction: a new Flex context vs the others), and
      every leave-one-reference-out subset; vs 639's mixed set. Gate sign accuracy on 636's h1 gate, paired vs CLEAN4,
      target bootstrap.
  (b) stability of I_A / I_B / I_C (refs = the other two A/B/C contexts [+ H1]) across disjoint NTC-guide halves of each
      context's controls: Spearman over genes, top-decile overlap, and the sign agreement of z(C4) + z(I) between halves.
Analysis only; < 4 GB (one context's control matrix at a time).
"""
from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402

RNG = np.random.default_rng(6220)
LAMS = (0.5, 1.0)


def abc_cpm(halves: bool = False):
    """Mean per-cell CPM of each context's controls (639's construction); with halves, also two disjoint NTC halves."""
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    genes = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
    out = {}
    for c in "ABC":
        a = ad.read_h5ad(ROOT / f"data/challenge/context_{c}.h5ad", backed="r")
        X = sp.csr_matrix(a.X[:])
        lib = np.asarray(X.sum(1)).ravel()
        N = sp.diags(1e6 / lib) @ X
        out[f"ctx{c}"] = np.asarray(N.mean(0)).ravel()
        if halves:
            ntc = a.obs["ntc_id"].astype(str).to_numpy()
            u = np.array(sorted(set(ntc)))
            RNG.shuffle(u)
            h1set = set(u[: len(u) // 2])
            m = np.array([x in h1set for x in ntc])
            out[f"ctx{c}_h1"] = np.asarray(N[m].mean(0)).ravel()
            out[f"ctx{c}_h2"] = np.asarray(N[~m].mean(0)).ravel()
        del X, N
    return genes, out


def zs(v):
    return (v - np.nanmean(v)) / np.nanstd(v)


def main() -> None:
    I = Inst()
    genes, abc = abc_cpm(halves=True)
    gpos = {g: i for i, g in enumerate(genes)}
    prof = {}
    for s, f in SURF.items():
        z = np.load(G480 / f, allow_pickle=True)
        prof[s] = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(float)))
    for k in ("ctxA", "ctxB", "ctxC"):
        prof[k] = dict(zip(genes, abc[k]))
    out = {"a_h1_reference_sets": [], "b_halves": {}}

    # ---- (a) h1 with different reference sets
    D = I.S["h1"]
    C4, _ = clean4_lolo(I, "h1")
    zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
    gate, dn = D["gate"], D["lfc"] < 0
    T = D["T"]

    def paired(score_g):
        """gate sign acc of score vs CLEAN4 on genes both cover; per-target sums for a bootstrap."""
        both = np.isfinite(zc) & np.isfinite(score_g)
        m = gate & both[None, :]
        okS = (((score_g > 0)[None, :] == dn) & m).sum(1).astype(float)
        okC = (((zc > 0)[None, :] == dn) & m).sum(1).astype(float)
        n = m.sum(1).astype(float)
        idx = np.flatnonzero(n > 0)
        bs = []
        for _ in range(2000):
            b = RNG.choice(idx, idx.size)
            bs.append((okS[b].sum() - okC[b].sum()) / n[b].sum())
        return dict(acc=float(okS.sum() / n.sum()), c4=float(okC.sum() / n.sum()), d=float((okS.sum() - okC.sum()) / n.sum()),
                    ci=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))], n=int(n.sum()))

    ref_sets = {"639 mixed (hct116, hek293t, k562, A, B, C)": ["hct116", "hek293t", "k562", "ctxA", "ctxB", "ctxC"],
                "Flex only {A, B, C} (board / day-of form)": ["ctxA", "ctxB", "ctxC"]}
    for drop in ("ctxA", "ctxB", "ctxC"):
        ref_sets[f"Flex, leave {drop} out"] = [k for k in ("ctxA", "ctxB", "ctxC") if k != drop]
    for name, refs in ref_sets.items():
        common = sorted(set(prof["h1"]).intersection(*(set(prof[r]) for r in refs)))
        Lh = np.log1p(np.array([prof["h1"][g] for g in common]))
        Lr = np.mean([np.log1p(np.array([prof[r][g] for g in common])) for r in refs], 0)
        Iv = I.gene_vec("h1", common, Lh - Lr)
        zi = zs(Iv)
        row = {"refs": name, "n_common": len(common)}
        for lam in LAMS:
            row[f"lam{lam}"] = paired(zc + lam * zi)
        row["identity_alone"] = paired(zi)
        out["a_h1_reference_sets"].append(row)
        print(f"[a] {name:44s} genes {len(common):5d} | " + "  ".join(
            f"lam{lam}: {row[f'lam{lam}']['d']:+.4f} [{row[f'lam{lam}']['ci'][0]:+.4f},{row[f'lam{lam}']['ci'][1]:+.4f}]" for lam in LAMS)
            + f"  | I alone {row['identity_alone']['acc']:.4f} vs C4 {row['identity_alone']['c4']:.4f}", flush=True)

    # ---- (b) stability of I for each A/B/C context across disjoint NTC halves (refs = other two contexts, full controls)
    c4g = {}
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(cp["bank_names"].astype(str))
    P = cp["P_banks"].astype(float)[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    c4map = dict(zip(cp["genes"].astype(str), P))
    c4v = np.array([c4map.get(g, np.nan) for g in genes])
    zc4 = (c4v - 0.5) / np.nanstd(c4v - 0.5)
    expressed = np.ones(len(genes), bool)
    for c in "ABC":
        others = [f"ctx{o}" for o in "ABC" if o != c]
        Lref = np.mean([np.log1p(abc[o]) for o in others], 0)
        I1 = np.log1p(abc[f"ctx{c}_h1"]) - Lref
        I2 = np.log1p(abc[f"ctx{c}_h2"]) - Lref
        Ifull = np.log1p(abc[f"ctx{c}"]) - Lref
        # restrict to genes the DE scorer would test in this context (> 5 CPM in controls)
        keep = abc[f"ctx{c}"] > 5
        from scipy.stats import spearmanr
        rho = float(spearmanr(I1[keep], I2[keep]).correlation)
        q1 = I1[keep] >= np.quantile(I1[keep], 0.9); q2 = I2[keep] >= np.quantile(I2[keep], 0.9)
        top_overlap = float((q1 & q2).sum() / q1.sum())
        res = {"genes_tested": int(keep.sum()), "spearman_I_halves": rho, "top_decile_overlap": top_overlap}
        for lam in LAMS:
            s1 = zc4 + lam * zs(I1); s2 = zc4 + lam * zs(I2)
            ok = keep & np.isfinite(s1) & np.isfinite(s2)
            res[f"sign_agree_lam{lam}"] = float(((s1[ok] > 0) == (s2[ok] > 0)).mean())
            s0 = zc4 + lam * zs(Ifull)
            ok0 = keep & np.isfinite(s0)
            res[f"frac_flipped_vs_C4_lam{lam}"] = float(((s0[ok0] > 0) != (zc4[ok0] > 0)).mean())
        out["b_halves"][c] = res
        print(f"[b] context {c}: genes {res['genes_tested']}, Spearman(I half1, half2) {rho:.4f}, top-decile overlap "
              f"{top_overlap:.3f}, score sign agreement λ0.5 {res['sign_agree_lam0.5']:.4f} λ1 {res['sign_agree_lam1.0']:.4f}; "
              f"genes flipped vs CLEAN4 λ0.5 {res['frac_flipped_vs_C4_lam0.5']:.3f} λ1 {res['frac_flipped_vs_C4_lam1.0']:.3f}", flush=True)
    (WORK / "stab622.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
