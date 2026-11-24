"""411 — Does CD4 (primary T cells) improve DE SIGNS for Jurkat (context A's line) beyond the three cancer donors?

Off-panel test on Nadig 2025's Jurkat essential screen (the challenge panel has no essential genes, so this is the only
same-line truth we have). Truth = Nadig Jurkat delta rows. Predictors, each line's rows mean-removed over the evaluated
targets and divided by that line's median row norm:
    C       mean of the available cancer donors (k562_gwps, orion_hct116, orion_hek293t; >= 2 required)
    C+λX    C + λ * X, X = CD4 (mean of available conditions) or the HepG2 control line (Nadig HepG2, same screen family)
Metrics per target, averaged: |t|-weighted sign accuracy on the truth's top-K |delta| genes (K = 50, 100, 200), and cosine.
Bootstrap over targets for the λ = 0 vs best-λ contrast. Writes work/results411.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent / "work"
OFFICIAL = [l.strip() for l in (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:] if l.strip()]
OFF = {g: i for i, g in enumerate(OFFICIAL)}
BANKS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
    "jurkat": ROOT / "data/donors/nadig_jurkat",
    "hepg2": ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2",
    "cd4_rest": WORK / "cd4_jurkat/cd4_rest",
    "cd4_stim8hr": WORK / "cd4_jurkat/cd4_stim8hr",
    "cd4_stim48hr": WORK / "cd4_jurkat/cd4_stim48hr",
}


def load(name, targets):
    d = BANKS[name]
    genes = [l.strip() for l in (d / "genes.txt").read_text().splitlines() if l.strip()]
    perts = [l.strip() for l in (d / "perts.txt").read_text().splitlines() if l.strip()]
    X = np.load(d / "delta.npy", mmap_mode="r")
    pi = {p: i for i, p in enumerate(perts)}
    src, dst = [], []
    for j, g in enumerate(genes):
        k = OFF.get(g)
        if k is not None:
            src.append(j); dst.append(k)
    out = np.zeros((len(targets), len(OFFICIAL)), dtype=np.float64)
    have = np.zeros(len(targets), dtype=bool)
    cov = np.zeros(len(OFFICIAL), dtype=bool); cov[dst] = True
    for r, t in enumerate(targets):
        i = pi.get(t)
        if i is not None:
            out[r, dst] = np.asarray(X[i], dtype=np.float64)[src]
            have[r] = True
    return out, have, cov


def prep(rows, have):
    """mean-remove over the carried targets, divide by the median row norm."""
    R = rows.copy()
    if have.sum() > 1:
        R[have] -= R[have].mean(0, keepdims=True)
        R[have] /= np.median(np.linalg.norm(R[have], axis=1))
    R[~have] = 0.0
    return R


def main():
    targets = [l.strip() for l in (BANKS["jurkat"] / "perts.txt").read_text().splitlines() if l.strip()]
    T, hT, covT = load("jurkat", targets)
    lines = {n: load(n, targets) for n in ("k562_gwps", "orion_hct116", "orion_hek293t", "hepg2", "cd4_rest", "cd4_stim8hr", "cd4_stim48hr")}
    cancer = ["k562_gwps", "orion_hct116", "orion_hek293t"]
    cd4c = ["cd4_rest", "cd4_stim8hr", "cd4_stim48hr"]
    # common gene axis: covered by the truth, every cancer donor, CD4 and HepG2; target's own gene dropped per row
    cov = covT.copy()
    for n in cancer + ["hepg2", "cd4_rest"]:
        cov &= lines[n][2]
    n_c = sum(lines[n][1].astype(int) for n in cancer)
    has_cd4 = np.any([lines[n][1] for n in cd4c], axis=0)
    has_h = lines["hepg2"][1]
    keep = hT & (n_c >= 2) & has_cd4 & has_h
    idx = np.flatnonzero(keep)
    print(f"targets: truth {hT.sum()}, >=2 cancer donors {(hT & (n_c >= 2)).sum()}, +CD4 {(hT & (n_c >= 2) & has_cd4).sum()}, "
          f"+HepG2 (evaluated) {keep.sum()}; genes on the common axis {cov.sum()}", flush=True)
    sub = [targets[i] for i in idx]
    Tm = prep(T[idx][:, cov], np.ones(len(idx), bool))
    P = {}
    for n in cancer + ["hepg2"]:
        rows, have, _ = lines[n]
        P[n] = (prep(rows[idx][:, cov], have[idx]), have[idx])
    cd4rows = np.zeros((len(idx), cov.sum())); cd4n = np.zeros(len(idx))
    for n in cd4c:
        rows, have, _ = lines[n]
        R = prep(rows[idx][:, cov], have[idx]); cd4rows += R; cd4n += have[idx]
    CD4 = cd4rows / np.maximum(cd4n, 1)[:, None]
    C = sum(P[n][0] for n in cancer) / np.maximum(sum(P[n][1].astype(int) for n in cancer), 1)[:, None]
    H = P["hepg2"][0]
    genes_cov = [OFFICIAL[k] for k in np.flatnonzero(cov)]
    gpos = {g: k for k, g in enumerate(genes_cov)}
    own = np.array([gpos.get(t, -1) for t in sub])

    def metrics(M, K):
        acc, cos = [], []
        for r in range(len(sub)):
            t = Tm[r].copy(); p = M[r].copy()
            if own[r] >= 0:
                t[own[r]] = 0.0; p[own[r]] = 0.0
            top = np.argsort(-np.abs(t))[:K]
            w = np.abs(t[top])
            acc.append(float((w * (np.sign(p[top]) == np.sign(t[top]))).sum() / w.sum()))
            cos.append(float(p @ t / (np.linalg.norm(p) * np.linalg.norm(t) + 1e-300)))
        return np.array(acc), np.array(cos)

    res = {"n_targets": len(sub), "n_genes": int(cov.sum()), "rows": {}}
    rng = np.random.default_rng(411)
    base = {K: metrics(C, K) for K in (50, 100, 200)}
    for xname, X in (("cd4", CD4), ("hepg2", H)):
        for lam in (0.0, 0.25, 0.5, 1.0, 2.0, 1e9):
            M = X.copy() if lam >= 1e9 else C + lam * X
            row = {}
            for K in (50, 100, 200):
                a, c = metrics(M, K)
                d = a - base[K][0]
                bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(1000)]
                row[f"K{K}"] = dict(acc=float(a.mean()), dacc=float(d.mean()), ci=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))],
                                    cos=float(c.mean()))
            res["rows"][f"{xname} lam {lam:g}"] = row
            r100 = row["K100"]
            print(f"{xname:5s} lam {lam:>6g}: signacc@100 {r100['acc']:.4f} (Δ {r100['dacc']:+.4f} [{r100['ci'][0]:+.4f},{r100['ci'][1]:+.4f}])  "
                  f"@50 Δ {row['K50']['dacc']:+.4f}  @200 Δ {row['K200']['dacc']:+.4f}  cos {r100['cos']:+.4f}", flush=True)
    (WORK / "results411.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
