#!/usr/bin/env python3
"""625b — which 16-module PARTITION recovers most of a lab programme? (before any read; may become Amendment 2)

Partitions (truth-free; defined on the challenge gene axis where possible, mapped to Broad-16 genes by symbol):
  c4xexpr   the preregistered one: CLEAN4 P_down octile x control-expression half
  coresp16  16 k-means clusters of genes' co-response across public CRISPRi perturbations (K562 GWPS: 3,000 random
            perturbations, per-gene robust z, PCA 20, cosine k-means), i.e. biological programmes
  coresp8x  8 co-response clusters x expression half
  c4xcoresp CLEAN4 sign (2) x 8 co-response clusters
For each, sim625's procedure on the Broad-16 (12 working lines; reserved 4 never read): 40 random 3-line "A/B/C" triples per
held-out line, module a_M = |l|-weighted net agreement of CLEAN4's sign on the triple's gate genes, flip CLEAN4 in modules
with a_M < -0.1, read held-out gate-sign gain vs CLEAN4. Genes outside a partition keep CLEAN4. Analysis only, < 6 GB.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/622-identity-track"))
import broad622 as B  # noqa: E402

RNG = np.random.default_rng(6251)
BAR = 0.10


def coresp_clusters(k: int, seed: int = 0) -> dict[str, int]:
    d = ROOT / "data/donors/k562_gwps"
    genes = [x.strip() for x in (d / "genes.txt").read_text().split("\n") if x.strip()]
    D = np.load(d / "delta.npy", mmap_mode="r")
    rng = np.random.default_rng(6252 + seed)
    rows = np.sort(rng.choice(D.shape[0], 3000, replace=False))
    X = np.asarray(D[rows], dtype=np.float64)                       # perts x genes
    med = np.median(X, 0); mad = np.median(np.abs(X - med), 0) * 1.4826
    Z = (X - med) / np.where(mad > 1e-9, mad, np.nan)
    Z = np.nan_to_num(np.clip(Z, -10, 10)).T                        # genes x perts
    Z -= Z.mean(1, keepdims=True)
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    E = U[:, :20] * S[:20]
    E /= np.maximum(np.linalg.norm(E, axis=1, keepdims=True), 1e-12)
    # cosine k-means (spherical), 20 restarts
    best, bestj = None, -np.inf
    for r in range(20):
        C = E[rng.choice(E.shape[0], k, replace=False)]
        for _ in range(50):
            lab = np.argmax(E @ C.T, 1)
            C = np.stack([E[lab == j].mean(0) if (lab == j).any() else E[rng.integers(E.shape[0])] for j in range(k)])
            C /= np.maximum(np.linalg.norm(C, axis=1, keepdims=True), 1e-12)
        j = float(np.sum(np.max(E @ C.T, 1)))
        if j > bestj:
            bestj, best = j, lab
    return dict(zip(genes, best.tolist()))


def main() -> None:
    gl, LFC, cl = B.load_matrix(B.P / "lfc_matrix.csv", 1)
    gf, FDR, cf = B.load_matrix(B.P / "fdr_matrix.csv", 1)
    gp, PB, cp = B.load_matrix(B.P / "pseudobulk_sum_by_ko.csv", 0)
    lines = sorted({ln for ln, _ in cl}); assert len(lines) == 12 and not (set(lines) & B.RESERVED)
    c4 = np.load(B.ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    nm = list(c4["bank_names"].astype(str))
    Pc = c4["P_banks"].astype(float)[[nm.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    cmap = dict(zip(c4["genes"].astype(str), Pc)); C4 = np.array([cmap.get(g, np.nan) for g in gl]); c4dn = C4 > 0.5
    all_t = set(k for _, k in cl)
    gpos = {g: i for i, g in enumerate(gp)}
    rowmap = np.array([gpos.get(g, -1) for g in gl])
    ctrl = {}
    for ln in lines:
        cc = [j for j, (l2, k) in enumerate(cp) if l2 == ln and k not in all_t]
        s = np.nansum(PB[:, cc], 1); v = 1e6 * s / np.nansum(s)
        ctrl[ln] = np.where(rowmap >= 0, v[np.maximum(rowmap, 0)], np.nan)
    cr16 = coresp_clusters(16); cr8 = coresp_clusters(8)
    m16 = np.array([cr16.get(g, -1) for g in gl]); m8 = np.array([cr8.get(g, -1) for g in gl])
    octile = np.digitize(C4, np.nanquantile(C4, np.linspace(0.125, 0.875, 7)))
    gated = (FDR < 0.05) & np.isfinite(LFC) & (LFC != 0)
    colline = np.array([ln for ln, _ in cl]); tgt = np.array([k for _, k in cl])
    notown = gl[:, None] != tgt[None, :]
    agree = np.where(gated & notown, np.where((LFC < 0) == c4dn[:, None], 1.0, -1.0), 0.0)
    absl = np.where(gated & notown, np.abs(LFC), 0.0)
    okc = np.isfinite(C4)
    names = ("c4xexpr", "coresp16", "coresp8x", "c4xcoresp")
    res = {n: [] for n in names}
    cover = {}
    for L in lines:
        others = [o for o in lines if o != L]
        jl = np.flatnonzero(colline == L)
        gL = gated[:, jl] & notown[:, jl]; tdnL = LFC[:, jl] < 0
        for _ in range(40):
            tri = list(RNG.choice(others, 3, replace=False))
            ex = np.nanmean([np.log1p(ctrl[o]) for o in tri], 0)
            half = (ex > np.nanmedian(ex)).astype(int)
            parts = {"c4xexpr": np.where(okc & np.isfinite(ex), octile * 2 + half, -1),
                     "coresp16": np.where(okc & (m16 >= 0), m16, -1),
                     "coresp8x": np.where(okc & (m8 >= 0) & np.isfinite(ex), m8 * 2 + half, -1),
                     "c4xcoresp": np.where(okc & (m8 >= 0), m8 * 2 + c4dn.astype(int), -1)}
            jt = np.flatnonzero(np.isin(colline, tri))
            num = (agree[:, jt] * absl[:, jt]).sum(1); den = absl[:, jt].sum(1)
            m = gL & okc[:, None]
            base = ((c4dn[:, None] == tdnL) & m).sum()
            for n in names:
                mod = parts[n]
                cover[n] = float((mod >= 0).mean())
                flip = np.zeros(16, bool)
                for k in range(16):
                    mk = mod == k
                    if den[mk].sum() > 0 and num[mk].sum() / den[mk].sum() < -BAR:
                        flip[k] = True
                fg = np.where(mod >= 0, flip[np.maximum(mod, 0)], False)
                sgn = np.where(fg, ~c4dn, c4dn)
                res[n].append(float((((sgn[:, None] == tdnL) & m).sum() - base) / max(m.sum(), 1)))
    out = {}
    for n in names:
        v = np.array(res[n]).reshape(12, 40)
        per = v.mean(1)
        out[n] = dict(mean=float(v.mean()), worst_line=float(per.min()), lines_pos=int((per > 0).sum()),
                      frac_triples_neg=float((v < 0).mean()), gene_coverage=cover[n])
        print(f"{n:10s} Δ gate-sign {out[n]['mean']:+.4f}  worst line {out[n]['worst_line']:+.4f}  lines>0 {out[n]['lines_pos']}/12  "
              f"triples<0 {out[n]['frac_triples_neg']:.2f}  gene coverage {cover[n]:.3f}", flush=True)
    (HERE / "work/partition625.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
