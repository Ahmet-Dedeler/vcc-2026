#!/usr/bin/env python3
"""633j — a second CLEAN cancer surface for the UP-dose question: Nadig 2025 Jurkat (CRISPRi, 10x; context A is Jurkat).

480's scorer-rule reference DE (gates480.py: Wilcoxon rank-sum vs the harness's own control cells, BH per perturbation
over genes > 5 CPM in control, lfc = log2 ratio of arithmetic mean CPM), adapted to this file's DENSE X (262,956 x 8,882
float32, old-style categorical obs). The panel's 300 targets are absent from this essential-gene screen (0/300), so the
surface is 300 random perturbations (seed 0) with >= 50 cells whose target K562 GWPS also perturbs (so the fold model's
donor consensus has a row). All non-targeting cells are the controls.
Writes work/gate-nadig_jurkat.npz in 480's format (targets, genes, ctrl_cpm, lib_med, lib_pert_med, ncell, n_ctrl, p, q,
lfc). Known answer: target 0's p / q / lfc are recomputed from its raw rows with scipy's Mann-Whitney (exact formula
differs only in the tie/continuity handling): see the KA block (reader bit-identity, >= 99 % call agreement, lfc).
Peak memory ~3-4 GB (the kept rows, float32).   uv run python experiments/633-updose-folds/gates633j.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import h5py
import numpy as np
from scipy.stats import false_discovery_control, mannwhitneyu, norm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/092-call-set-truncation"))
from de_kernel import ctrl_tiesums, wilcoxon_z  # noqa: E402

F = ROOT / "data/external/nadig2025/jurkat_raw_singlecell.h5ad"
N_T, MIN_CELLS, SEED, BLK, CHUNK = 300, 50, 0, 8192, 1000


def main() -> int:
    t0 = time.time()
    h = h5py.File(F, "r")
    lab = h["obs"]["__categories"]["gene"][:].astype(str)[h["obs"]["gene"][:]]
    vg = h["var"]["gene_name"]
    genes = (h["var"]["__categories"]["gene_name"][:].astype(str)[vg[:]] if "gene_name" in h["var"]["__categories"]
             else vg[:].astype(str))
    X = h["X"]
    n, G = X.shape
    gw = [l.strip() for l in open(ROOT / "data/donors/k562_gwps/perts.txt")]
    u, c = np.unique(lab, return_counts=True)
    cand = sorted(p for p, k in zip(u, c) if p != "non-targeting" and k >= MIN_CELLS and p in set(gw) and p in set(genes))
    rng = np.random.default_rng(SEED)
    perts = sorted(rng.choice(cand, size=min(N_T, len(cand)), replace=False).tolist())
    cidx = np.flatnonzero(lab == "non-targeting")
    rows = {p: np.flatnonzero(lab == p) for p in perts}
    keep = np.sort(np.concatenate([cidx] + list(rows.values())))
    print(f"[jurkat] {n} cells x {G} genes; candidates {len(cand)}, chose {len(perts)}; controls {cidx.size}; "
          f"kept rows {keep.size}", flush=True)
    M = np.empty((keep.size, G), dtype=np.float32)
    k = 0
    for lo in range(0, n, BLK):
        hi = min(lo + BLK, n)
        sel = keep[(keep >= lo) & (keep < hi)]
        if sel.size:
            M[k:k + sel.size] = X[lo:hi][sel - lo]
            k += sel.size
    h.close()
    assert k == keep.size
    if not np.all(M == np.floor(M)) or M.min() < 0:
        raise SystemExit("X is not non-negative integer counts")
    pos = np.full(n, -1)
    pos[keep] = np.arange(keep.size)
    lib = M.sum(1).astype(np.float64)
    lib[lib == 0] = 1.0
    c_pos = pos[cidx]
    w = np.zeros(keep.size)
    w[c_pos] = 1e6 / lib[c_pos] / c_pos.size
    cpm_all = M.T.astype(np.float64) @ w
    tested = np.flatnonzero(cpm_all > 5.0)
    print(f"[jurkat] loaded ({time.time() - t0:.0f}s); {tested.size} genes > 5 CPM; median ctrl lib "
          f"{np.median(lib[c_pos]):.0f}", flush=True)
    Gt = tested.size
    P = np.ones((len(perts), Gt))
    LFC = np.zeros((len(perts), Gt), dtype=np.float32)
    for lo in range(0, Gt, CHUNK):
        cols = tested[lo:lo + CHUNK]
        D = M[:, cols].astype(np.float64) * (1e6 / lib)[:, None]
        Cc = D[c_pos]
        cmean = Cc.mean(0)
        cs = np.ascontiguousarray(np.sort(np.log1p(Cc).T, axis=1))
        tie = ctrl_tiesums(cs)
        for i, p in enumerate(perts):
            Bv = D[pos[rows[p]]]
            z = wilcoxon_z(cs, tie, np.ascontiguousarray(np.log1p(Bv).T))
            P[i, lo:lo + cols.size] = 2 * norm.sf(np.abs(z))
            LFC[i, lo:lo + cols.size] = np.log2((Bv.mean(0) + 1e-9) / (cmean + 1e-9))
        print(f"  genes {lo + cols.size}/{Gt} ({time.time() - t0:.0f}s)", flush=True)
    Q = np.vstack([false_discovery_control(P[i], method="bh") for i in range(len(perts))])
    gk = genes[tested]
    # KA (independent of the block reader): target 0's rows and 2,000 control rows re-read by index from the file must
    # equal the loaded matrix bit for bit; scipy's Mann-Whitney on them must agree with the kernel's BH calls on >= 99 %
    # of genes (the kernel is 480/092's scanpy-identical one; tie handling may differ at the margin); lfc exact.
    p0 = perts[0]
    rs = np.sort(np.concatenate([rows[p0], np.random.default_rng(1).choice(cidx, 2000, replace=False)]))
    with h5py.File(F, "r") as h2:
        R2 = h2["X"][rs]
    reader_ok = np.array_equal(R2, M[pos[rs]])
    Dr = M[:, tested].astype(np.float64) * (1e6 / lib)[:, None]
    Cr, Br = np.log1p(Dr[c_pos]), np.log1p(Dr[pos[rows[p0]]])
    pk = mannwhitneyu(Br, Cr, alternative="two-sided", method="asymptotic", use_continuity=False, axis=0).pvalue
    qk = false_discovery_control(np.nan_to_num(pk, nan=1.0), method="bh")
    lk = np.log2((Dr[pos[rows[p0]]].mean(0) + 1e-9) / (Dr[c_pos].mean(0) + 1e-9))
    agree = float(((qk < 0.05) == (Q[0] < 0.05)).mean())
    lfc_ok = np.allclose(lk, LFC[0], atol=1e-5)
    print(f"[jurkat] KA target {p0}: reader bit-identical {reader_ok}; BH call agreement with scipy {agree:.4f} "
          f"({int((Q[0] < 0.05).sum())} kernel calls, {int((qk < 0.05).sum())} scipy); lfc equal {lfc_ok}", flush=True)
    if not (reader_ok and agree >= 0.99 and lfc_ok):
        raise SystemExit("KA FAIL")
    ncell = np.array([rows[p].size for p in perts])
    sig = [(Q[i] < 0.05) & np.isfinite(LFC[i]) & (LFC[i] != 0) & (gk != p) for i, p in enumerate(perts)]
    sz = np.array([s.sum() for s in sig])
    lib_pert = np.median(lib[pos[np.concatenate(list(rows.values()))]])
    print(f"[jurkat] gate sizes median {np.median(sz):.0f} mean {sz.mean():.0f}; >=10: {(sz >= 10).sum()}/{len(sz)}; "
          f"median lib ctrl {np.median(lib[c_pos]):.0f} pert {lib_pert:.0f}", flush=True)
    np.savez_compressed(HERE / "work/gate-nadig_jurkat.npz", targets=np.array(perts), genes=gk, ctrl_cpm=cpm_all[tested],
                        lib_med=np.median(lib[c_pos]), lib_pert_med=lib_pert, ncell=ncell, n_ctrl=c_pos.size,
                        p=P.astype(np.float32), q=Q.astype(np.float32), lfc=LFC)
    print(f"[jurkat] wrote work/gate-nadig_jurkat.npz ({time.time() - t0:.0f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
