#!/usr/bin/env python3
"""517 — 462's rule applied to the hook's source: is Arc's H1 programme's accuracy on h1-test's genes inflated by the
CONTROL cells the programme and the truth share? (h1_2025_full's programme is taken against all 38,176 public H1
controls, which include h1-test's 18,400; the harness truth is taken against those 18,400.)

Every public control cell gets a half (md5 of its cell id, so a cell has the same half in the raw files and in h1-test):
  programme  P_X = P_down over the 215 non-h1-test perturbations of lp(pert sums) - lp(control half X sums)   (X = A, B, all)
  truth      T_X = log2 arithmetic-mean-CPM fold change of each h1-test target vs h1-test's controls in half X
  shared   = P_A vs T_A (and P_B vs T_B);   disjoint = P_A vs T_B (and P_B vs T_A);   full = P_all vs T_all (the harness)
Gate = the scorer's reference table (q < .05, finite lfc != 0, own gene out). Regions: 462f's block (2,200 top control-CPM
kept genes, panel targets excluded) and the rest (the hook's territory), plus fid-like (all block genes, gate or not).
KAs: P_all's signs == the h1_2025_full bank's programme (h1-test rows out); T_all == the reference lfc (float32 accumulation).
Streams everything (< 2 GB).  usage: uv run python experiments/517-flex-programme/halves517.py
"""
from __future__ import annotations

import csv
import glob
import hashlib
import json
import time
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FILES = [ROOT / "data/external/vcc2025/train/adata_Training.h5ad", ROOT / "data/external/vcc2025/validation/adata_Validation.h5ad",
         ROOT / "data/external/vcc2025/test/adata_Test.h5ad"]
H1T = ROOT / "data/harness/h1-test.h5ad"
TS, CHUNK, EPS = 5e4, 4000, 1e-9
t0 = time.time()
log = lambda *x: print("[517h]", *x, f"({time.time() - t0:.0f}s)", flush=True)   # noqa: E731


def cat(n):
    return n["categories"][:].astype(str)[n["codes"][:]] if isinstance(n, h5py.Group) else n[:].astype(str)


def ids_of(f):
    ix = f["obs/_index"]
    return (ix[:] if not isinstance(ix, h5py.Group) else ix["values"][:]).astype(str)


def half_of(ids):
    return np.array([int(hashlib.md5(s.encode()).hexdigest()[:8], 16) & 1 for s in ids], dtype=np.int8)


def buffers(f, ipd):
    """Preallocated read buffers for the largest CHUNK-row slice (read_direct): per-chunk allocations made the first
    version of this script peak at 6.7 GB RSS on macOS (freed chunks were not reused); buffered, it stays < 1 GB."""
    mx = int(max(ipd[min(s + CHUNK, ipd.size - 1)] - ipd[s] for s in range(0, ipd.size - 1, CHUNK)))
    return np.empty(mx, f["X/data"].dtype), np.empty(mx, f["X/indices"].dtype)


def main():
    z = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=False)
    genes, keep, targets = z["genes"].astype(str), np.asarray(z["keep"], bool), list(z["targets"].astype(str))
    q_ref, lfc_ref, bctrl = np.asarray(z["q"], float), np.asarray(z["lfc"], float), np.asarray(z["b_ctrl"], float)
    z.close()
    kg = genes[keep]; kcol = np.flatnonzero(keep)
    h1set = set(targets)
    # ---------------------------------------------------------------- 1. programme sums from the raw public files
    sums: dict[str, np.ndarray] = {}
    ctrl = np.zeros((3, genes.size))            # all, half A, half B
    seen: set[str] = set()
    for path in FILES:
        with h5py.File(path, "r") as f:
            gv = f["var/_index"]
            gv = (gv[:] if not isinstance(gv, h5py.Group) else gv["values"][:]).astype(str)
            assert np.array_equal(gv, genes), f"{path.name}: gene axis differs from the harness"
            tv, ids = cat(f["obs/target_gene"]), ids_of(f)
            is_c = tv == "non-targeting"
            dup = is_c & np.isin(ids, list(seen)) if seen else np.zeros(tv.size, bool)
            seen |= set(ids[is_c])
            hv = np.where(is_c, half_of(np.where(is_c, ids, "")), -1)
            use_p = ~is_c & ~np.isin(tv, list(h1set))
            ipd = f["X/indptr"][:].astype(np.int64)
            dbuf, ibuf = buffers(f, ipd)
            for s in range(0, tv.size, CHUNK):
                e = min(tv.size, s + CHUNK)
                a, b = int(ipd[s]), int(ipd[e])
                f["X/data"].read_direct(dbuf, np.s_[a:b], np.s_[0:b - a]); f["X/indices"].read_direct(ibuf, np.s_[a:b], np.s_[0:b - a])
                M = sp.csr_matrix((dbuf[:b - a], ibuf[:b - a], ipd[s:e + 1] - a), shape=(e - s, genes.size), copy=False)
                sl = slice(s, e)
                cm = is_c[sl] & ~dup[sl]
                if cm.any():
                    for h in (0, 1):
                        sel = cm & (hv[sl] == h)
                        if sel.any():
                            v = np.asarray(M[np.flatnonzero(sel)].sum(0)).ravel()
                            ctrl[0] += v; ctrl[1 + h] += v
                pm = use_p[sl]
                if pm.any():
                    loc = np.flatnonzero(pm)
                    for t in np.unique(tv[sl][loc]):
                        rows = loc[tv[sl][loc] == t]
                        sums.setdefault(t, np.zeros(genes.size))
                        sums[t] += np.asarray(M[rows].sum(0)).ravel()
            log(f"{path.name}: {tv.size} cells; {int(is_c.sum())} controls ({int(dup.sum())} already counted)")
    perts = sorted(sums)
    log(f"programme rows {len(perts)} (non-h1-test), controls {len(seen)} (A {int(ctrl[1].sum() > 0)})")
    lp = lambda c: np.log1p(TS * c / c.sum())   # noqa: E731

    def pdown(cX):
        dn = np.zeros(genes.size); nz = np.zeros(genes.size); lc = lp(cX)
        for p in perts:
            d = (lp(sums[p]) - lc).astype(np.float32)     # the bank's frame (delta.npy is float32)
            dn += d < 0; nz += d != 0
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)[keep]
    PF, PA, PB = pdown(ctrl[0]), pdown(ctrl[1]), pdown(ctrl[2])
    # KA 1: P_all's signs == the h1_2025_full bank's programme with h1-test's rows out (467's arc_pdown)
    d = ROOT / "data/donors/h1_2025_full"
    delta = np.load(d / "delta.npy", mmap_mode="r"); bp = [l.strip() for l in open(d / "perts.txt")]
    rows = [i for i, p in enumerate(bp) if p not in h1set]
    Xb = np.asarray(delta[rows], dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        Pbank = np.where((Xb != 0).sum(0) > 0, (Xb < 0).sum(0) / np.maximum((Xb != 0).sum(0), 1), np.nan)[keep]
    del Xb
    agree = float(np.mean(np.sign(0.5 - PF[np.isfinite(PF)]) == np.sign(0.5 - Pbank[np.isfinite(PF)])))
    maxd = float(np.nanmax(np.abs(PF - Pbank)))
    log(f"KA1: programme(all controls) vs the bank's: max |dP| {maxd:.2e}, sign agreement {agree:.6f} ({len(rows)} bank rows)")
    assert len(rows) == len(perts) and maxd < 0.01, "KA1 failed"
    # ---------------------------------------------------------------- 2. truth from h1-test vs its control halves
    with h5py.File(H1T, "r") as f:
        tv, ids = cat(f["obs/target_gene"]), ids_of(f)
        is_c = tv == "non-targeting"
        hv = np.where(is_c, half_of(np.where(is_c, ids, "")), -1)
        tidx = {t: i for i, t in enumerate(targets)}
        S_t = np.zeros((len(targets), kg.size)); n_t = np.zeros(len(targets))
        S_c = np.zeros((3, kg.size)); n_c = np.zeros(3)
        ipd = f["X/indptr"][:].astype(np.int64)
        dbuf, ibuf = buffers(f, ipd)
        for s in range(0, tv.size, CHUNK):
            e = min(tv.size, s + CHUNK)
            a, b = int(ipd[s]), int(ipd[e])
            f["X/data"].read_direct(dbuf, np.s_[a:b], np.s_[0:b - a]); f["X/indices"].read_direct(ibuf, np.s_[a:b], np.s_[0:b - a])
            M = sp.csr_matrix((dbuf[:b - a].astype(np.float64), ibuf[:b - a], ipd[s:e + 1] - a), shape=(e - s, genes.size))
            lib = np.asarray(M.sum(1)).ravel()
            C = (sp.diags(1e6 / lib) @ M)[:, kcol].tocsr()
            sl = slice(s, e)
            for h in (-1, 0, 1):
                sel = is_c[sl] & ((hv[sl] == h) if h >= 0 else True)
                if sel.any():
                    S_c[h + 1] += np.asarray(C[np.flatnonzero(sel)].sum(0)).ravel(); n_c[h + 1] += sel.sum()
            pm = np.flatnonzero(~is_c[sl])
            for t in np.unique(tv[sl][pm]):
                rr = pm[tv[sl][pm] == t]
                S_t[tidx[t]] += np.asarray(C[rr].sum(0)).ravel(); n_t[tidx[t]] += rr.size
    Mt = S_t / n_t[:, None]
    mF, mA, mB = (S_c[i] / n_c[i] for i in range(3))
    lg = lambda m, c: np.log2((m + EPS) / (c + EPS))   # noqa: E731
    TF, TA, TB = lg(Mt, mF), lg(Mt, mA), lg(Mt, mB)
    fin = np.isfinite(lfc_ref) & np.isfinite(TF)
    ka2 = float(np.nanmax(np.abs(TF[fin] - lfc_ref[fin])))
    log(f"KA2: truth(all h1-test controls) vs the reference lfc: max |d| {ka2:.2e}, median {np.median(np.abs(TF[fin] - lfc_ref[fin])):.1e}; "
        f"control halves {int(n_c[1])} / {int(n_c[2])}")
    assert ka2 < 1e-3, "KA2 failed"
    # ---------------------------------------------------------------- 3. accuracies
    e_ = np.expm1(bctrl[keep]); cpm = 1e6 * e_ / e_.sum()
    order = np.argsort(-cpm)
    block = np.zeros(kg.size, bool); block[[i for i in order if kg[i] not in h1set][:2200]] = True
    own = np.array([[g == t for g in kg] for t in targets])
    gate = (q_ref < 0.05) & np.isfinite(lfc_ref) & (lfc_ref != 0) & ~own
    c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
    offi = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
    pos = {g: i for i, g in enumerate(offi)}
    c4p = np.asarray(c4["P_down"], float)
    PC4 = np.array([c4p[pos[g]] if g in pos else np.nan for g in kg])

    def acc(P, T, region, gated=True, tau=0.0):
        s = np.where(np.isfinite(P) & (np.abs(P - 0.5) >= tau), np.where(P > 0.5, -1.0, 1.0), np.nan)
        hit = n = 0
        for t in range(len(targets)):
            m = region & np.isfinite(s) & np.isfinite(T[t]) & (T[t] != 0) & ~own[t]
            m &= gate[t] if gated else True
            hit += int((np.sign(T[t][m]) == s[m]).sum()); n += int(m.sum())
        return round(hit / max(n, 1), 4), n
    res = {"ka1_max_dP": maxd, "ka1_sign_agree": agree, "ka2_max_dlfc": ka2}
    for name, (region, gated, tau) in {"gate_block": (block, True, 0.0), "gate_nonblock": (~block, True, 0.0),
                                       "gate_nonblock_tau05": (~block, True, 0.05), "fidlike_block": (block, False, 0.0),
                                       "all_nonblock_tau05": (~block, False, 0.05)}.items():
        row = {"full (P_all vs T_all)": acc(PF, TF, region, gated, tau),
               "shared (P_A vs T_A)": acc(PA, TA, region, gated, tau), "shared (P_B vs T_B)": acc(PB, TB, region, gated, tau),
               "disjoint (P_A vs T_B)": acc(PA, TB, region, gated, tau), "disjoint (P_B vs T_A)": acc(PB, TA, region, gated, tau),
               "CLEAN4 vs T_A": acc(PC4, TA, region, gated, tau), "CLEAN4 vs T_B": acc(PC4, TB, region, gated, tau)}
        res[name] = row
        sh = np.mean([row["shared (P_A vs T_A)"][0], row["shared (P_B vs T_B)"][0]])
        dj = np.mean([row["disjoint (P_A vs T_B)"][0], row["disjoint (P_B vs T_A)"][0]])
        c4m = np.mean([row["CLEAN4 vs T_A"][0], row["CLEAN4 vs T_B"][0]])
        log(f"{name:22s} full {row['full (P_all vs T_all)'][0]:.4f} | shared {sh:.4f} disjoint {dj:.4f} (inflation {sh - dj:+.4f}) "
            f"| CLEAN4 {c4m:.4f} | n {row['full (P_all vs T_all)'][1]}")
        res[name + "_summary"] = {"shared": sh, "disjoint": dj, "inflation": sh - dj, "clean4": c4m}
    (HERE / "work/halves517.json").write_text(json.dumps(res, indent=1))
    log("done")


if __name__ == "__main__":
    main()
