"""295 item 3 — THE SINGLE ACCESS to reserve R3.  One pass per file, everything computed here.

Preregistration `experiments/295-reserve-r3-claim/PREREGISTRATION.md`
sha256 a8a49c637c6c24e9503513d84576030c13b911b0dbc8bb89d455751ff78b8a7d, appended to
docs/CONFIRMATION-LEDGER.md as claim C5 at 2026-09-06T16:09:10Z, BEFORE this script ran.

Per file, in one streaming pass pair (pass 1: per-cell library size; pass 2: the axis columns):
  * the residual response of every gene with >= 50 cells that is in the held map
  * its split-half reliability (8 draws)
  * its DE count (Welch t on log1p(CPM) vs the file's controls, BH q<0.05, control CPM > 5)
  * cells per gene, median UMIs per cell per gene and per file
  * KA-5 (computed library size == obs.ncounts) and KA-6 (controls-vs-controls split ~ 0)

Writes work/access295.npz (per-row arrays), work/access295.json (per-file summary + KA rows).
"""
from __future__ import annotations

import csv
import json
import zlib
import resource
import sys
import time
from pathlib import Path

import h5py
import numpy as np
from scipy import stats as sstats

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/295-reserve-r3-claim"
WORK = HERE / "work"

TS = 5e4
SEED = 20260906
NDRAW = 8
ROWBLK = 8192
CTRL_CPM_MIN = 5.0
FDR_Q = 0.05

PREREG_SHA = "a8a49c637c6c24e9503513d84576030c13b911b0dbc8bb89d455751ff78b8a7d"


def check_prereg() -> None:
    import hashlib
    h = hashlib.sha256((HERE / "PREREGISTRATION.md").read_bytes()).hexdigest()
    if h != PREREG_SHA:
        raise SystemExit(f"preregistration digest changed ({h}) — abort")
    print(f"preregistration digest {h} OK", flush=True)


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30


def _decode(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a)
    if a.dtype.kind in ("S", "O"):
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in a])
    return a.astype(str)


def read_str_col(f: h5py.File, path: str) -> np.ndarray:
    d = f[path]
    if isinstance(d, h5py.Group):
        cats = _decode(d["categories"])
        codes = np.asarray(d["codes"], np.int64)
        # code -1 = pandas NaN.  278's obs scan serialised it as the EMPTY string (verified on
        # Adamson: 2,613 cells, key ""), which 292's frozen CTRL_PAT matches as a control.
        # Never index `cats` with -1.
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "")
    return _decode(d)


def bh(p: np.ndarray, q: float) -> np.ndarray:
    """Benjamini-Hochberg rejection mask."""
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    out = np.zeros(p.shape, bool)
    pv = p[ok]
    if pv.size == 0:
        return out
    order = np.argsort(pv)
    ps = pv[order]
    m = ps.size
    thresh = q * np.arange(1, m + 1) / m
    below = ps <= thresh
    k = np.flatnonzero(below)
    if k.size:
        cut = k[-1]
        rej = np.zeros(m, bool)
        rej[order[:cut + 1]] = True
        out[ok] = rej
    return out


# ------------------------------------------------------------------ the two streaming passes
def pass1_totals(f: h5py.File, enc: str, n_cells: int) -> np.ndarray:
    """Per-cell library size, summed over ALL genes."""
    X = f["X"]
    data = X["data"]
    if enc == "csr_matrix":
        indptr = np.asarray(X["indptr"], np.int64)
        tot = np.zeros(n_cells, np.float64)
        for lo in range(0, n_cells, ROWBLK):
            hi = min(lo + ROWBLK, n_cells)
            a, b = int(indptr[lo]), int(indptr[hi])
            if b <= a:
                continue
            d = np.asarray(data[a:b], np.float64)
            counts = np.diff(indptr[lo:hi + 1])
            rid = np.repeat(np.arange(hi - lo, dtype=np.int64), counts)
            tot[lo:hi] = np.bincount(rid, weights=d, minlength=hi - lo)
        return tot
    # csc: need the cell index of every nonzero
    idxs = X["indices"]
    tot = np.zeros(n_cells, np.float64)
    nnz = int(data.shape[0])
    step = 1 << 24
    for lo in range(0, nnz, step):
        hi = min(lo + step, nnz)
        d = np.asarray(data[lo:hi], np.float64)
        i = np.asarray(idxs[lo:hi], np.int64)
        tot += np.bincount(i, weights=d, minlength=n_cells)
    return tot


def pass2_accumulate(f: h5py.File, enc: str, n_cells: int, cols: np.ndarray,
                     acc_of_col: np.ndarray, nA: int,
                     col2ax: np.ndarray, grp: np.ndarray, half: np.ndarray,
                     tot: np.ndarray, n_grp: int):
    """Accumulate, over the axis columns only:
         Sraw[grp, a]           sum of raw counts
         Slog[grp, a], Qlog     sum and sum of squares of log1p(1e6 * v / tot)
         Half[grp, draw, h, a]  sum of raw counts per split-half
    """
    Sraw = np.zeros((n_grp, nA), np.float64)
    Slog = np.zeros((n_grp, nA), np.float64)
    Qlog = np.zeros((n_grp, nA), np.float64)
    Half = np.zeros((n_grp, NDRAW, 2, nA), np.float64)
    X = f["X"]
    data, idxs = X["data"], X["indices"]
    indptr = np.asarray(X["indptr"])
    inv_tot = np.where(tot > 0, 1.0 / np.maximum(tot, 1.0), 0.0)

    if enc == "csc_matrix":
        for k, c in enumerate(cols):
            j = int(acc_of_col[k])
            a, b = int(indptr[c]), int(indptr[c + 1])
            if b <= a:
                continue
            v = np.asarray(data[a:b], np.float64)
            cell = np.asarray(idxs[a:b], np.int64)
            g = grp[cell]
            keep = g >= 0
            if not keep.any():
                continue
            v, cell, g = v[keep], cell[keep], g[keep]
            lg = np.log1p(1e6 * v * inv_tot[cell])
            Sraw[:, j] += np.bincount(g, weights=v, minlength=n_grp)
            Slog[:, j] += np.bincount(g, weights=lg, minlength=n_grp)
            Qlog[:, j] += np.bincount(g, weights=lg * lg, minlength=n_grp)
            for d in range(NDRAW):
                key = g * 2 + half[cell, d]
                Half[:, d, :, j] += np.bincount(key, weights=v,
                                                minlength=2 * n_grp).reshape(n_grp, 2)
        return Sraw, Slog, Qlog, Half

    # csr
    for lo in range(0, n_cells, ROWBLK):
        hi = min(lo + ROWBLK, n_cells)
        a, b = int(indptr[lo]), int(indptr[hi])
        if b <= a:
            continue
        v = np.asarray(data[a:b], np.float64)
        cj = np.asarray(idxs[a:b], np.int64)
        counts = np.diff(indptr[lo:hi + 1])
        cell = np.repeat(np.arange(lo, hi, dtype=np.int64), counts)
        ax = col2ax[cj]
        g = grp[cell]
        keep = (ax >= 0) & (g >= 0)
        if not keep.any():
            continue
        v, cell, ax, g = v[keep], cell[keep], ax[keep], g[keep]
        lg = np.log1p(1e6 * v * inv_tot[cell])
        flat = g * nA + ax
        Sraw += np.bincount(flat, weights=v, minlength=n_grp * nA).reshape(n_grp, nA)
        Slog += np.bincount(flat, weights=lg, minlength=n_grp * nA).reshape(n_grp, nA)
        Qlog += np.bincount(flat, weights=lg * lg, minlength=n_grp * nA).reshape(n_grp, nA)
        for d in range(NDRAW):
            key = (g * 2 + half[cell, d]) * nA + ax
            Half[:, d, :, :] += np.bincount(key, weights=v,
                                            minlength=n_grp * 2 * nA).reshape(n_grp, 2, nA)
    return Sraw, Slog, Qlog, Half


def prof(S: np.ndarray, T: float) -> np.ndarray:
    return np.log1p(TS * S / T) if T > 0 else np.zeros_like(S)


def cos(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else np.nan


def main() -> None:
    check_prereg()
    t0 = time.time()
    meta = json.loads((WORK / "meta295.json").read_text())
    mp = np.load(WORK / "map295.npz", allow_pickle=True)
    U, avail = mp["U"], mp["avail"]                      # (28, 174, 2739), (28, 174)
    mgenes = [str(g) for g in mp["genes"]]
    mkeys = [str(k) for k in mp["keys"]]
    axis = [str(g) for g in mp["axis"]]
    axis_pos = {g: i for i, g in enumerate(axis)}
    gpos = {g: i for i, g in enumerate(mgenes)}
    rng_master = np.random.default_rng(SEED)

    rows, per_file, ka5, ka6 = [], {}, {}, {}
    for stem in sorted(meta["files"]):
        M = meta["files"][stem]
        t1 = time.time()
        p = ROOT / M["path"]
        sel = M["genes"]
        lab_of = M["labels_per_gene"]
        ctrl_labels = set(M["control_labels"])
        excl_banks = set(M["excluded_banks"])
        bank_ok = np.asarray([k not in excl_banks for k in mkeys])

        with h5py.File(p, "r") as f:
            enc = M["structure"]["X"]["attrs"]["encoding-type"]
            n_cells, n_var = M["structure"]["X"]["attrs"]["shape"]
            lab = read_str_col(f, "obs/perturbation")
            lab = np.asarray([s.strip() for s in lab])
            ncounts = np.asarray(f["obs/ncounts"], np.float64)
            vidx = M["structure"]["var_attrs"]["_index"]
            vsym = read_str_col(f, f"var/{vidx}")

            # ---- groups: 0 = control, 1..K = the selected genes, -1 = unused ------------------
            grp = np.full(n_cells, -1, np.int64)
            grp[np.isin(lab, list(ctrl_labels))] = 0
            for gi, g in enumerate(sel, start=1):
                grp[np.isin(lab, lab_of[g])] = gi
            n_grp = len(sel) + 1
            sizes = np.bincount(grp[grp >= 0], minlength=n_grp)
            assert int(sizes[0]) == M["n_control_cells"], (stem, sizes[0], M["n_control_cells"])
            for gi, g in enumerate(sel, start=1):
                assert int(sizes[gi]) == M["cells_per_gene"][g], (stem, g)

            # ---- the column set: the map axis genes this file carries -------------------------
            present = {s for s in vsym if s in axis_pos}
            names = [g for g in axis if g in present]
            name_index = {g: i for i, g in enumerate(names)}
            cols = np.asarray([i for i, s in enumerate(vsym) if s in present], np.int64)
            acc_of_col = np.asarray([name_index[vsym[i]] for i in cols], np.int64)
            ax_of_name = np.asarray([axis_pos[g] for g in names], np.int64)
            col2ax = np.full(n_var, -1, np.int64)
            col2ax[cols] = acc_of_col            # duplicate symbols sum into one accumulator
            assert len(vsym) == n_var, (stem, len(vsym), n_var)

            # ---- half assignment: per cell, per draw ------------------------------------------
            rng = np.random.default_rng(zlib.crc32(stem.encode()) ^ SEED)
            half = np.zeros((n_cells, NDRAW), np.int64)
            for gi in range(n_grp):
                idx = np.flatnonzero(grp == gi)
                if idx.size == 0:
                    continue
                for d in range(NDRAW):
                    perm = rng.permutation(idx.size)
                    h = np.zeros(idx.size, np.int64)
                    h[perm[idx.size // 2:]] = 1
                    half[idx, d] = h

            print(f"[{stem}] {enc} {n_cells}x{n_var}; {len(sel)} genes; axis genes {len(names)} "
                  f"in {len(cols)} columns",
                  flush=True)
            tot = pass1_totals(f, enc, n_cells)
            d5 = float(np.max(np.abs(tot - ncounts)))
            rel5 = float(np.max(np.abs(tot - ncounts) / np.maximum(ncounts, 1.0)))
            ka5[stem] = dict(max_abs_diff=d5, max_rel_diff=rel5, pass_=bool(rel5 < 1e-6))
            print(f"    KA-5 library size vs obs.ncounts: max|diff| {d5:.3g}, "
                  f"max rel {rel5:.3g} -> {'PASS' if ka5[stem]['pass_'] else 'FAIL'}", flush=True)
            if not ka5[stem]["pass_"]:
                raise SystemExit(f"KA-5 FAILED on {stem} — abort")

            Sraw, Slog, Qlog, Half = pass2_accumulate(
                f, enc, n_cells, cols, acc_of_col, len(names), col2ax, grp, half, tot, n_grp)

        # ---- group library totals ------------------------------------------------------------
        T = np.zeros(n_grp)
        umed = np.zeros(n_grp)
        for gi in range(n_grp):
            m = grp == gi
            T[gi] = float(tot[m].sum())
            umed[gi] = float(np.median(tot[m])) if m.any() else np.nan
        file_umi = float(np.median(tot[grp >= 0]))

        # ---- responses ------------------------------------------------------------------------
        P = np.vstack([prof(Sraw[gi], T[gi]) for gi in range(n_grp)])       # (n_grp, nA)
        D = P - P[0][None, :]                                              # vs the file's control
        R = D[1:] - D[1:].mean(axis=0, keepdims=True)                      # residual, genes only

        # split halves (each half centred on its own across-gene mean)
        rsplit = np.full(len(sel), np.nan)
        Th = np.zeros((n_grp, NDRAW, 2))
        for gi in range(n_grp):
            m = np.flatnonzero(grp == gi)
            for d in range(NDRAW):
                for h in (0, 1):
                    Th[gi, d, h] = float(tot[m[half[m, d] == h]].sum())
        for d in range(NDRAW):
            for h in (0, 1):
                Ph = np.vstack([prof(Half[gi, d, h], Th[gi, d, h]) for gi in range(n_grp)])
                Dh = Ph - Ph[0][None, :]
                Rh = Dh[1:] - Dh[1:].mean(axis=0, keepdims=True)
                if h == 0:
                    R0 = Rh
                else:
                    c = np.asarray([cos(R0[i], Rh[i]) for i in range(len(sel))])
                    rsplit = c if d == 0 else rsplit + c
        rsplit = rsplit / NDRAW

        # ---- KA-6: controls vs controls --------------------------------------------------------
        Pc0 = prof(Half[0, 0, 0], Th[0, 0, 0])
        Pc1 = prof(Half[0, 0, 1], Th[0, 0, 1])
        Dnull = Pc0 - Pc1

        names_g = names
        # ---- DE counts (Welch t on log1p CPM, BH q<0.05, control pseudobulk CPM > 5) -----------
        n_of = sizes.astype(float)
        mu = Slog / np.maximum(n_of[:, None], 1)
        var = (Qlog - n_of[:, None] * mu ** 2) / np.maximum(n_of[:, None] - 1, 1)
        var = np.maximum(var, 0.0)
        cpm_ctrl = 1e6 * Sraw[0] / max(T[0], 1.0)
        gate = cpm_ctrl > CTRL_CPM_MIN
        de_count = np.zeros(len(sel), int)
        for gi in range(1, n_grp):
            n1, n2 = n_of[gi], n_of[0]
            se2 = var[gi] / n1 + var[0] / n2
            with np.errstate(invalid="ignore", divide="ignore"):
                t = (mu[gi] - mu[0]) / np.sqrt(se2)
                df = se2 ** 2 / ((var[gi] / n1) ** 2 / max(n1 - 1, 1)
                                 + (var[0] / n2) ** 2 / max(n2 - 1, 1))
            pv = np.full(len(names_g), np.nan)
            ok = gate & np.isfinite(t) & (se2 > 0) & np.isfinite(df) & (df > 0)
            pv[ok] = 2 * sstats.t.sf(np.abs(t[ok]), df[ok])
            de_count[gi - 1] = int(bh(pv, FDR_Q).sum())

        # ---- the outcome ------------------------------------------------------------------------
        ka6_cos = []
        for i, g in enumerate(sel):
            j = gpos[g]
            use = bank_ok & avail[:, j]
            nb = int(use.sum())
            if nb == 0:
                rows.append(dict(file=stem, gene=g, A=np.nan, A_onT=np.nan, n_banks=0,
                                 cells=int(sizes[i + 1]), guides=M["guides_per_gene"][g],
                                 r_split=float(rsplit[i]), de=int(de_count[i]),
                                 umi_gene=float(umed[i + 1]), umi_file=file_umi,
                                 line=M["line"], modality=M["modality"]))
                continue
            C = U[use][:, j, :].mean(axis=0)                 # consensus on the FULL 2,739 axis
            c_sub = C[ax_of_name]                            # this file's columns
            keep = np.ones(len(names), bool)
            if g in axis_pos:
                w = np.flatnonzero(ax_of_name == axis_pos[g])
                keep[w] = False
            a_on = cos(R[i], c_sub)
            a_off = cos(R[i][keep], c_sub[keep])
            ka6_cos.append(cos(Dnull[keep], c_sub[keep]))
            rows.append(dict(file=stem, gene=g, A=a_off, A_onT=a_on, n_banks=nb,
                             cells=int(sizes[i + 1]), guides=M["guides_per_gene"][g],
                             r_split=float(rsplit[i]), de=int(de_count[i]),
                             umi_gene=float(umed[i + 1]), umi_file=file_umi,
                             line=M["line"], modality=M["modality"]))
        med_null = float(np.nanmedian(ka6_cos)) if ka6_cos else np.nan
        ka6[stem] = dict(median_cos=med_null, n=len(ka6_cos),
                         p10=float(np.nanpercentile(ka6_cos, 10)) if ka6_cos else np.nan,
                         p90=float(np.nanpercentile(ka6_cos, 90)) if ka6_cos else np.nan,
                         pass_=bool(abs(med_null) < 0.05))
        A_here = np.asarray([r["A"] for r in rows if r["file"] == stem], float)
        per_file[stem] = dict(
            n_genes=len(sel), enc=enc, n_axis_cols=len(names), file_umi_median=file_umi,
            control_cells=int(sizes[0]), median_cells=float(np.median(sizes[1:])),
            median_r_split=float(np.nanmedian(rsplit)), median_de=float(np.median(de_count)),
            median_A=float(np.nanmedian(A_here)), ka5=ka5[stem], ka6=ka6[stem],
            seconds=time.time() - t1)
        print(f"    KA-6 controls-vs-controls median cos vs consensus {med_null:+.4f} "
              f"-> {'PASS' if ka6[stem]['pass_'] else 'FAIL'}", flush=True)
        print(f"    median cells {np.median(sizes[1:]):.0f}, UMI/cell {file_umi:.0f}, "
              f"r_split {np.nanmedian(rsplit):+.4f}, DE {np.median(de_count):.0f}, "
              f"median A {np.nanmedian(A_here):+.4f}   ({time.time()-t1:.0f}s, "
              f"peak {peak_gib():.2f} GiB)", flush=True)
        if not ka6[stem]["pass_"]:
            raise SystemExit(f"KA-6 FAILED on {stem} — abort")
        del Sraw, Slog, Qlog, Half, P, D, R, tot, grp, half

    # ---- save --------------------------------------------------------------------------------
    keysf = ["file", "gene", "A", "A_onT", "n_banks", "cells", "guides", "r_split", "de",
             "umi_gene", "umi_file", "line", "modality"]
    with open(WORK / "access295.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keysf)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    out = dict(n_rows=len(rows), per_file=per_file, ka5=ka5, ka6=ka6,
               seconds=time.time() - t0, peak_rss_gib=peak_gib(),
               prereg_sha256=PREREG_SHA)
    (WORK / "access295.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"\n{len(rows)} rows written; {time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB",
          flush=True)


if __name__ == "__main__":
    main()
