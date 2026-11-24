"""516 step 5: denoise each donor row before combining (directions 2 and 4 of the mandate).

For each bank on the surface's gene axis (K562-GWPS, Orion HCT116, Orion HEK293T, CD4 per condition; never the surface's
own line; the surface's targets are never in the basis):
  A     = 6,000 random other perturbations' log2 ratios (CD4: DESeq2 log2FC), centred by the per-gene median (the bank
          programme); the surface's targets are never among them
  SVD   randomized, top k = 400 (Halko, 2 power iterations), noise power per component from the residual energy
  wiener(x) = med + V diag(s^2 / (s^2 + n^2)) V^T xc + rho (xc - V V^T xc)       (rho = residual weight, 515's shape)
  knn(x)    = med + xc + lam * mean_{q in N_k(p)} (A_q)   N_k = the k most correlated other perturbations in this bank
            (related-target borrowing: complex / pathway partners' rows shrink a noisy row toward its module)
The denoised rows replace the raw rows in the Stouffer (same per-gene mad and per-row tau as the raw row, so only the
sign structure changes) and are priced on H1 (vs production) and on the folds (vs the production analog).
    PYTHONDONTWRITEBYTECODE=1 uv run python denoise516.py build <surface>     # h1 | orion_hct116 | orion_hek293t | k562_gwps
    PYTHONDONTWRITEBYTECODE=1 uv run python denoise516.py price
"""
from __future__ import annotations

import json
import sys
import time

import h5py
import numpy as np

from common516 import CD4_H5AD, ROOT, WORK, known_answers, load_h1, price, summary_line
from sources516 import TS, surface_axis

K_SVD, OVERS, NPI = 400, 30, 2
N_BASIS = 6000          # random other perturbations per bank (basis and kNN pool); keeps every process < 2 GB
RHOS = (0.0, 0.5)
KNN = ((10, 0.5), (30, 1.0))
BANKS = ["k562_gwps", "orion_hct116", "orion_hek293t"]
CD4C = ["Rest", "Stim8hr", "Stim48hr"]


def rsvd(A, k, seed=516):
    rng = np.random.default_rng(seed)
    Om = rng.standard_normal((A.shape[1], k + OVERS)).astype(np.float32)
    Y = A @ Om
    for _ in range(NPI):
        Q, _ = np.linalg.qr(Y)
        Y = A @ (Q.T @ A).T                           # (Q^T A)^T avoids a transposed copy of A
    Q, _ = np.linalg.qr(Y)
    B = Q.T @ A
    Ub, s, Vt = np.linalg.svd(B, full_matrices=False)
    tot = sum(float(np.dot(A[i:i + 1000].ravel(), A[i:i + 1000].ravel())) for i in range(0, A.shape[0], 1000))
    top = float((s[:k].astype(np.float64) ** 2).sum())
    rest_dims = max(min(A.shape) - k, 1)
    n2 = (tot - top) / rest_dims
    return s[:k], Vt[:k], n2


def denoise_rows(A, med, X, has):
    """A: other perts (centred, float32, genes on axis cols); X: surface rows (T x g_axis, raw, NaN ok)."""
    s, Vt, n2 = rsvd(A, K_SVD)
    w = (s ** 2 / (s ** 2 + n2)).astype(np.float32)
    Xc = np.nan_to_num(X - med[None, :]).astype(np.float32)
    P = Xc @ Vt.T
    proj = P @ Vt
    out = {}
    for rho in RHOS:
        out[f"wie{rho:g}"] = med[None, :] + (P * w[None, :]) @ Vt + rho * (Xc - proj)
    # kNN in the bank: correlation of each surface row with every other perturbation (centred rows), in row chunks
    Xn = Xc - Xc.mean(1, keepdims=True); Xn /= np.maximum(np.linalg.norm(Xn, axis=1, keepdims=True), 1e-9)
    C = np.empty((len(X), A.shape[0]), np.float32)
    for a in range(0, A.shape[0], 2000):
        An = A[a:a + 2000] - A[a:a + 2000].mean(1, keepdims=True)
        An /= np.maximum(np.linalg.norm(An, axis=1, keepdims=True), 1e-9)
        C[:, a:a + 2000] = Xn @ An.T
    for k, lam in KNN:
        idx = np.argpartition(-C, k, axis=1)[:, :k]
        nb = np.stack([A[idx[i]].mean(0) for i in range(len(X))])
        out[f"knn{k}_{lam:g}"] = med[None, :] + Xc + lam * nb
    info = dict(n_other=int(A.shape[0]), s_top=[float(x) for x in s[:5]], n2=float(n2),
                w_mean=float(w.mean()), w_first=float(w[0]), w_last=float(w[-1]),
                knn_r_top1_median=float(np.median(np.sort(C, 1)[:, -1])))
    for kk in out:
        out[kk] = np.where(has[:, None], out[kk], np.nan).astype(np.float32)
    return out, info


class NpyRows:
    """Read whole rows of a C-order 2-D .npy with seek/readinto (no mmap, so RSS stays at the working set)."""
    def __init__(self, path):
        import numpy.lib.format as fmt
        self.f = open(path, "rb")
        v = fmt.read_magic(self.f)
        shape, fortran, dtype = (fmt.read_array_header_1_0 if v == (1, 0) else fmt.read_array_header_2_0)(self.f)
        assert not fortran and len(shape) == 2
        self.off, self.shape, self.dtype = self.f.tell(), shape, np.dtype(dtype)
        self.rowb = shape[1] * self.dtype.itemsize

    def rows(self, idx):
        out = np.empty((len(idx), self.shape[1]), self.dtype)
        for j, i in enumerate(idx):
            self.f.seek(self.off + int(i) * self.rowb)
            self.f.readinto(memoryview(out[j]).cast("B"))
        return out


def bank_matrix(name, genes, exclude):
    base = ROOT / "data/donors" / name
    P = [l.strip() for l in open(base / "perts.txt")]
    Gb = [l.strip() for l in open(base / "genes.txt")]
    gcol = {}
    for i, g in enumerate(Gb):
        if g and g not in gcol:
            gcol[g] = i
    b = np.asarray(np.load(base / "b_ctrl.npy"), float)
    cols = np.array([gcol.get(g, -1) for g in genes]); okc = cols >= 0; cc = cols[okc]
    rows = np.array([i for i, p in enumerate(P) if p not in exclude])
    rows = np.sort(np.random.default_rng(516).choice(rows, min(N_BASIS, len(rows)), replace=False))
    A = np.empty((len(rows), okc.sum()), np.float32)
    bc = b[cc][None, :]
    ctrl = np.expm1(bc) / TS * 1e6
    rd = NpyRows(base / "delta.npy")                  # plain reads, no mmap residency
    for a in range(0, len(rows), 500):
        r = rows[a:a + 500]
        d = rd.rows(r)[:, cc].astype(np.float64)
        pert = np.maximum(np.expm1(bc + d), 0) / TS * 1e6
        A[a:a + len(r)] = np.log2((pert + 1.0) / (ctrl + 1.0))
    return A, okc


def cd4_matrix(cond, genes, exclude):
    f = h5py.File(CD4_H5AD, "r")
    obs = f["obs"]
    o = obs["target_contrast_gene_name"]; tname = o["categories"][:].astype(str)[o["codes"][:]]
    o = obs["culture_condition"]; cnd = o["categories"][:].astype(str)[o["codes"][:]]
    gname = f["var"]["gene_name"][:].astype(str)
    gcol = {}
    for i, g in enumerate(gname):
        gcol.setdefault(g, i)
    cols = np.array([gcol.get(g, -1) for g in genes]); okc = cols >= 0; cc = cols[okc]
    rows = np.flatnonzero((cnd == cond) & ~np.isin(tname, list(exclude)))
    rows = np.sort(np.random.default_rng(516).choice(rows, min(N_BASIS, len(rows)), replace=False))
    A = np.empty((len(rows), okc.sum()), np.float32)
    for a in range(0, len(rows), 300):
        r = rows[a:a + 300]
        A[a:a + len(r)] = np.nan_to_num(f["layers"]["log_fc"][r, :][:, cc])
    f.close()
    return A, okc


def build(surf):
    t0 = time.time()
    targets, genes = surface_axis(surf)
    exclude = set(targets)
    z = np.load(WORK / f"src-{surf}.npz", allow_pickle=True)
    out, infos = {}, {}
    own = None if surf == "h1" else surf
    for bn in BANKS + [f"cd4_{c}" for c in CD4C]:
        if bn == own:
            continue
        if bn.startswith("cd4_"):
            A, okc = cd4_matrix(bn[4:], genes, exclude)
            X = z[f"{bn}__log_fc"].astype(np.float32)
        else:
            A, okc = bank_matrix(bn, genes, exclude)
            X = z[f"{bn}__x"].astype(np.float32)
        med = np.concatenate([np.median(A[:, j:j + 1000], 0) for j in range(0, A.shape[1], 1000)])
        A -= med[None, :]
        has = np.isfinite(X).any(1)
        res, info = denoise_rows(A, med, X[:, okc], has)
        del A
        for k, v in res.items():
            full = np.full((len(targets), len(genes)), np.nan, np.float32)
            full[:, okc] = v
            full[~np.isfinite(X)] = np.nan                  # keep the raw row's coverage exactly
            out[f"{bn}__{k}"] = full
        infos[bn] = info
        print(f"{surf}: {bn} basis {info['n_other']} rows, s[:3] {np.round(info['s_top'][:3], 1)}, noise/comp "
              f"{info['n2']:.2f}, wiener w first/last {info['w_first']:.2f}/{info['w_last']:.2f}, kNN top-1 r median "
              f"{info['knn_r_top1_median']:.3f} ({time.time() - t0:.0f}s)", flush=True)
    np.savez_compressed(WORK / f"den-{surf}.npz", **out)
    json.dump(infos, open(WORK / f"den-{surf}.json", "w"), indent=1)


def std_from(z, name, x, cap=3.0):
    return np.clip(x / (z[f"{name}__mad"][None, :] * z[f"{name}__tau"][:, None]), -cap, cap)


def stouffer_rows(z, rows):
    num = 0; den = 0
    for nm, x in rows:
        s = std_from(z, nm, x.astype(float))
        num = num + np.nan_to_num(s); den = den + np.isfinite(s)
    return np.where(den > 0, num / np.sqrt(np.maximum(den, 1)), np.nan)


def price_all():
    import fold516 as FD
    d = load_h1()
    got, ok, _ = known_answers(d)
    print(f"KNOWN ANSWER {'PASS' if ok else 'FAIL'}: fid {got['fid']!r} jac {got['jac']!r} nmae {got['nmae']!r} "
          f"reach {got['reach']!r}", flush=True)
    assert ok
    z = np.load(WORK / "src-h1.npz", allow_pickle=True); dn = np.load(WORK / "den-h1.npz")
    ess = np.zeros(len(d["targets"]), bool)
    for bn in ("k562_essential", "rpe1", "nadig_jurkat"):
        ess |= z[f"{bn}__has"]
    SUBS = {"all85": np.ones(len(d["targets"]), bool), "no_ess": ~ess}
    srcs = BANKS + [f"cd4_{c}" for c in CD4C]
    blk = np.zeros(len(d["genes"]), bool); blk[np.argsort(-d["ctrl"])[:2200]] = True
    variants = ["raw"] + [f"wie{r:g}" for r in RHOS] + [f"knn{k}_{l:g}" for k, l in KNN]
    rows = []

    def rowsof(zz, dd, names, v):
        return [(n, (zz[f"{n}__log_fc"] if n.startswith("cd4_") else zz[f"{n}__x"]) if v == "raw" else dd[f"{n}__{v}"])
                for n in names]
    # single-bank sign accuracy (K562 alone, as in 515) and the six-source Stouffer, per variant
    for v in variants:
        for label, names in (("K562 alone", ["k562_gwps"]), ("stz6", srcs)):
            Z = stouffer_rows(z, rowsof(z, dn, names, v))
            for thr, region, rm in ((0.0, "uncalled", ~d["calls"]), (0.5, "uncalled", ~d["calls"]),
                                    (0.0, "all", np.ones_like(d["calls"]))):
                C = np.where(rm & (np.abs(np.nan_to_num(Z)) >= thr), np.sign(np.nan_to_num(Z)), 0.0)
                r = price(d, C, f"H1 {label} [{v}] |Z|>={thr} {region}", subsets=SUBS, verbose=False)
                for sn in SUBS:
                    print(summary_line(r[sn]), flush=True)
                    rows.append(dict(surface="h1", label=label, variant=v, thr=thr, region=region, subset=sn,
                                     acc=r[sn]["acc"], members=r[sn]["members"]))
    for name, spec in FD.FOLDS.items():
        try:
            dnf = np.load(WORK / f"den-{name}.npz")
        except FileNotFoundError:
            print(f"(no den-{name}.npz; fold skipped)")
            continue
        Fo = FD.load_fold(name)
        donors = [x for x in spec["prod"] if x != spec["sib"]]
        A, blkf = FD.analogs(Fo, spec, donors)
        names = donors + [f"cd4_{c}" for c in CD4C]
        for v in variants:
            Z = stouffer_rows(Fo["z"], rowsof(Fo["z"], dnf, names, v))
            for thr in (0.0, 0.5):
                C = np.where(~blkf[None, :] & (np.abs(np.nan_to_num(Z)) >= thr), np.sign(np.nan_to_num(Z)), 0.0)
                FD.price_fold(Fo, A["base_mcd4+c4q80"], C, f"stz [{v}] |Z|>={thr} rest vs mcd4+c4q80", rows)
                FD.price_fold(Fo, A["base_mean"], C, f"stz [{v}] |Z|>={thr} rest vs mean", rows)
    json.dump(rows, open(WORK / "denoise516.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    if sys.argv[1] == "build":
        build(sys.argv[2])
    else:
        price_all()
