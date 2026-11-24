"""516: cache same-target source rows and per-bank noise models on a surface's (targets x genes) axis.

For every donor bank: the same-target row x (log2 pseudobulk ratio, pc 1 CPM; common480 conversion) and a rank-1
robust noise model built from the bank's OTHER perturbations (a random 2,000, never the surface's targets):
    med_g   per-gene median of x over perturbations  (the bank's generic programme at gene g)
    mad_g   per-gene median |x - med_g|               (gene noise scale: most perturbations do not move most genes)
    tau_p   per-row median |x_pg - med_g| / mad_g     (row noise scale: few cells -> noisy row)
so z_pg = (x_pg - med_g) / (mad_g * tau_p) is a bank-internal standardised deviation and x_pg / (mad_g * tau_p) its
uncentred analogue. CD4 (GWCD4i, DESeq2): per condition log_fc, zscore (Wald), plus the same med/mad of log_fc over a
random 2,000 targets. Banks h1_2025 / h1_2025_full are never loaded for the H1 surface (oracle leakage).

    PYTHONDONTWRITEBYTECODE=1 uv run python sources516.py h1
    PYTHONDONTWRITEBYTECODE=1 uv run python sources516.py orion_hct116 | orion_hek293t | k562_gwps
"""
from __future__ import annotations

import sys
import time

import h5py
import numpy as np

from common516 import CD4_H5AD, ROOT, WORK, load_h1

TS = 5e4
BANKS = ["k562_gwps", "orion_hct116", "orion_hek293t", "viperturb_k562", "k562_essential", "rpe1", "nadig_jurkat"]
import os
N_SAMPLE, SEED = int(os.environ.get("SRC_N_SAMPLE", 2000)), 516


def l2_rows(M, b, rows, cols, pc=1.0):
    """log2((pert CPM + pc) / (ctrl CPM + pc)) for bank rows `rows` at bank columns `cols` (float32)."""
    d = np.asarray(M[rows][:, cols], dtype=np.float64)
    bc = b[cols][None, :]
    pert = np.maximum(np.expm1(bc + d), 0) / TS * 1e6
    ctrl = np.expm1(bc) / TS * 1e6
    return np.log2((pert + pc) / (ctrl + pc)).astype(np.float32)


def bank(name, targets, genes, exclude):
    base = ROOT / "data/donors" / name
    P = [l.strip() for l in open(base / "perts.txt")]
    Gb = [l.strip() for l in open(base / "genes.txt")]
    pi = {p: i for i, p in enumerate(P)}
    gcol = {}
    for i, g in enumerate(Gb):
        if g and g not in gcol:
            gcol[g] = i
    M = np.load(base / "delta.npy", mmap_mode="r")
    b = np.asarray(np.load(base / "b_ctrl.npy"), float)
    cols = np.array([gcol.get(g, -1) for g in genes]); okc = cols >= 0
    cc = cols[okc]
    rng = np.random.default_rng(SEED)
    pool = [i for i, p in enumerate(P) if p not in exclude]
    samp = np.sort(rng.choice(pool, min(N_SAMPLE, len(pool)), replace=False))
    S = np.concatenate([l2_rows(M, b, samp[i:i + 250], cc) for i in range(0, len(samp), 250)])
    med = np.median(S, 0)
    S -= med[None, :]
    np.abs(S, out=S)
    mad = np.median(S, 0)
    mad = np.maximum(mad, np.quantile(mad[mad > 0], 0.01))
    del S
    X = np.full((len(targets), len(genes)), np.nan, np.float32)
    tau = np.full(len(targets), np.nan)
    has = np.array([t in pi for t in targets])
    if has.any():
        rr = np.array([pi[t] for t in targets[has]])
        order = np.argsort(rr)
        R = l2_rows(M, b, rr[order], cc)
        Rin = np.empty_like(R); Rin[order] = R
        Xh = np.full((has.sum(), len(genes)), np.nan, np.float32); Xh[:, okc] = Rin
        X[has] = Xh
        tau[has] = np.median(np.abs(Rin - med[None, :]) / mad[None, :], 1)
    MED = np.full(len(genes), np.nan, np.float32); MED[okc] = med
    MAD = np.full(len(genes), np.nan, np.float32); MAD[okc] = mad
    CPM = np.full(len(genes), np.nan, np.float32); CPM[okc] = np.expm1(b[cc]) / TS * 1e6
    return dict(x=X, med=MED, mad=MAD, tau=tau.astype(np.float32), cpm=CPM, has=has, n_bank=len(P))


def cd4(targets, genes, exclude):
    f = h5py.File(CD4_H5AD, "r")
    obs = f["obs"]

    def cat(k):
        o = obs[k]
        return o["categories"][:].astype(str)[o["codes"][:]]
    tname = cat("target_contrast_gene_name"); cond = cat("culture_condition")
    gname = f["var"]["gene_name"][:].astype(str)
    gcol = {}
    for i, g in enumerate(gname):
        gcol.setdefault(g, i)
    cols = np.array([gcol.get(g, -1) for g in genes]); okc = cols >= 0
    out = {}
    rng = np.random.default_rng(SEED)
    for c in sorted(set(cond)):
        rows = {k: int(np.flatnonzero((tname == t) & (cond == c))[0]) for k, t in enumerate(targets)
                if ((tname == t) & (cond == c)).any()}
        cand = np.flatnonzero((cond == c) & ~np.isin(tname, list(exclude)))
        samp = np.sort(rng.choice(cand, min(N_SAMPLE, len(cand)), replace=False))
        S = np.concatenate([np.nan_to_num(f["layers"]["log_fc"][samp[i:i + 250], :][:, cols[okc]].astype(np.float32))
                            for i in range(0, len(samp), 250)])
        med = np.median(S, 0); S -= med[None, :]; np.abs(S, out=S); mad = np.median(S, 0)
        mad = np.maximum(mad, np.quantile(mad[mad > 0], 0.01))
        del S
        res = {}
        ks = sorted(rows, key=lambda k: rows[k]); ridx = [rows[k] for k in ks]
        for L in ("log_fc", "zscore", "lfcSE", "baseMean"):
            A = np.full((len(targets), len(genes)), np.nan, np.float32)
            if ridx:
                A[np.array(ks)[:, None], np.flatnonzero(okc)[None, :]] = f["layers"][L][ridx, :][:, cols[okc]]
            res[L] = A
        tau = np.full(len(targets), np.nan)
        if ks:
            tau[np.array(ks)] = np.nanmedian(np.abs(res["log_fc"][np.array(ks)][:, okc] - med[None, :]) / mad[None, :], 1)
        res["tau"] = tau.astype(np.float32)
        MED = np.full(len(genes), np.nan, np.float32); MED[okc] = med
        MAD = np.full(len(genes), np.nan, np.float32); MAD[okc] = mad
        res["med"], res["mad"] = MED, MAD
        for mk in ("n_cells_target", "donor_correlation_all_mean", "guide_correlation_all", "n_total_de_genes"):
            v = np.full(len(targets), np.nan)
            if ridx:
                v[np.array(ks)] = obs[mk][ridx].astype(float)
            res["meta_" + mk] = v
        out[c] = res
    f.close()
    return out


def surface_axis(name):
    if name == "h1":
        d = load_h1()
        return d["targets"], d["genes"]
    z = np.load(ROOT / f"experiments/480-nmae-instrument/work/gate-{name}.npz", allow_pickle=True)
    return z["targets"].astype(str), z["genes"].astype(str)


def main():
    surf = sys.argv[1]
    t0 = time.time()
    targets, genes = surface_axis(surf)
    exclude = set(targets)
    out = {"targets": targets, "genes": genes}
    for bn in BANKS:
        r = bank(bn, targets, genes, exclude)
        for k, v in r.items():
            out[f"{bn}__{k}"] = v
        print(f"{surf}: {bn} rows {r['has'].sum()}/{len(targets)}, genes on axis {np.isfinite(r['mad']).sum()}, "
              f"tau median {np.nanmedian(r['tau']):.3f} ({time.time() - t0:.0f}s)", flush=True)
    c = cd4(targets, genes, exclude)
    for cn, r in c.items():
        for k, v in r.items():
            out[f"cd4_{cn}__{k}"] = v
        print(f"{surf}: cd4 {cn} rows {np.isfinite(r['log_fc']).any(1).sum()}/{len(targets)} "
              f"tau median {np.nanmedian(r['tau']):.3f} ({time.time() - t0:.0f}s)", flush=True)
    np.savez_compressed(WORK / f"src-{surf}.npz", **out)
    print(f"wrote work/src-{surf}.npz ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
