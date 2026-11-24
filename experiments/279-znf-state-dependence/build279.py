"""279 — the HELD-gene conservation map, using 276's frozen statistic unchanged.

KA-1 first (276's panel map reproduced from `conservation276.csv`), then the map over the
17,106 non-panel genes measured in >= 2 cell lines, its per-pair cosine matrix, the per-gene
response magnitude, the CCLE expression level and the target-shuffled null.

    PYTHONDONTWRITEBYTECODE=1 uv run python build279.py

Writes work/held279.npz, work/ka279.json, work/build279.json.  Read-only outside work/.
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/279-znf-state-dependence"
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/276-conservation-map"))
import lib276 as L  # noqa: E402
import inv279 as I  # noqa: E402

SEED = 20260906
N_NULL = 50
BLK = 8192


# ---------------------------------------------------------------- KA-1: 276's panel map
def ka1() -> dict:
    keys = list(L.SOURCES)
    axis = L.common_axis(keys)
    cross = L.pair_list(keys, cross_line_only=True)
    C = np.full((len(cross), len(L.PANEL)), np.nan)
    R, A = [], []
    for k in keys:
        r, a = L.load_source(k, axis)
        R.append(r); A.append(a)
    for p, (i, j) in enumerate(cross):
        both = A[i] & A[j]
        C[p, both] = L.cosine_rows(R[i][both].astype(np.float64), R[j][both].astype(np.float64))
    with np.errstate(invalid="ignore"):
        c_raw = np.nanmean(C, axis=0)
    ref = {}
    import csv as _csv
    with open(ROOT / "experiments/276-conservation-map/work/conservation276.csv") as fh:
        for row in _csv.DictReader(fh):
            ref[row["target"]] = float(row["C_raw"])
    got = np.asarray([c_raw[L.PANEL_IDX[t]] for t in L.PANEL])
    want = np.asarray([ref[t] for t in L.PANEL])
    d = float(np.nanmax(np.abs(got - want)))
    out = dict(n_axis=len(axis), n_cross_pairs=len(cross),
               median=float(np.nanmedian(c_raw)), share_le0=float(np.nanmean(c_raw <= 0)),
               max_abs_diff_vs_276=d, pass_=bool(d < 1e-6 and len(axis) == 2994
                                                 and len(cross) == 315))
    print(f"KA-1  axis {len(axis)} (want 2994), cross pairs {len(cross)} (want 315), "
          f"median {np.nanmedian(c_raw):+.7f} (want +0.0198672), "
          f"share<=0 {np.nanmean(c_raw <= 0):.3f} (want 0.160), "
          f"max|diff| vs conservation276.csv {d:.2e}  -> {'PASS' if out['pass_'] else 'FAIL'}")
    return out


# ---------------------------------------------------------------- the held map
def held_genes(per_src: dict) -> list[str]:
    panel = set(L.PANEL)
    official = set(L.OFFICIAL)
    lines: dict[str, set[str]] = {}
    for k, d in per_src.items():
        for g in d["perts"]:
            if g in panel or g not in official:
                continue
            lines.setdefault(g, set()).add(d["line"])
    return sorted(g for g, ls in lines.items() if len(ls) >= 2)


def residual_rows(rel: str, axis: list[str], genes: list[str], row_blk: int = 1024):
    """`lib276.load_source_path` in row blocks, with the delta memmap released afterwards.

    Identical arithmetic (KA-0 asserts it byte-for-byte); the block loop only keeps the peak RSS
    down — `D[rows][:, cols]` materialises the whole bank for the two 18k-gene Orion tables."""
    perts, bank_genes, p = L.read_bank(rel)
    gidx = {g: i for i, g in enumerate(bank_genes)}
    cols = np.asarray([gidx[g] for g in axis], dtype=np.int64)
    pidx = {t: i for i, t in enumerate(perts)}
    have = [(t, pidx[t]) for t in genes if t in pidx]
    if not have:
        return np.zeros(0, np.int64), np.zeros((0, len(axis)), np.float64)
    tidx = {t: i for i, t in enumerate(genes)}
    order = np.argsort([r for _, r in have])
    rows_sorted = np.asarray([have[i][1] for i in order], dtype=np.int64)
    idx = np.asarray([tidx[have[i][0]] for i in order], dtype=np.int64)
    dest = np.empty(idx.size, np.int64)             # write straight into gene-sorted order, so the
    dest[np.argsort(idx)] = np.arange(idx.size)     # 1 GiB block is never copied a second time
    D = np.load(p / "delta.npy", mmap_mode="r")
    blk = np.empty((rows_sorted.size, len(axis)), dtype=np.float64)
    for lo in range(0, rows_sorted.size, row_blk):
        r = rows_sorted[lo:lo + row_blk]
        blk[dest[lo:lo + r.size]] = np.asarray(D[r.min():r.max() + 1][:, cols],
                                               dtype=np.float64)[r - r.min()]
    D._mmap.close(); del D
    blk -= blk.mean(axis=0, keepdims=True)          # centre on this source's across-gene mean
    return np.sort(idx), blk


def load_unit(rel: str, axis: list[str], genes: list[str]):
    """residual rows (DEFINITION-276 s2) restricted to `genes`, unit-normalised, plus their index."""
    idx, blk = residual_rows(rel, axis, genes)
    if idx.size == 0:
        return idx, np.zeros((0, len(axis)), np.float32), np.zeros(0)
    rows = blk.astype(np.float32)
    del blk
    nrm = np.linalg.norm(rows, axis=1)
    ok = nrm > 0
    idx, rows, nrm = idx[ok], rows[ok], nrm[ok]
    rows /= nrm[:, None]
    return idx, rows, nrm


def ka0(axis: list[str], genes: list[str]) -> dict:
    """KA-0: the block reader reproduces lib276.load_source_path exactly on a mid-size bank."""
    rel = I.CANDIDATES["gse295214_h1"][0]
    R, A = L.load_source_path(rel, axis, targets=genes)
    idx, blk = residual_rows(rel, axis, genes)
    d = float(np.max(np.abs(R[idx].astype(np.float64) - blk))) if idx.size else np.nan
    same = bool(np.array_equal(np.flatnonzero(A), idx))
    print(f"KA-0  block reader vs lib276.load_source_path on gse295214_h1: "
          f"same index {same}, max|diff| {d:.2e} -> {'PASS' if same and d < 1e-6 else 'FAIL'}")
    return dict(same_index=same, max_abs_diff=d, pass_=bool(same and d < 1e-6))


def pair_cos(idx_i, U_i, idx_j, U_j):
    """cosine for the genes both sources carry; returns (gene indices, cosines)."""
    common = np.intersect1d(idx_i, idx_j, assume_unique=True)
    if common.size == 0:
        return common, np.zeros(0, np.float32)
    pi = np.searchsorted(idx_i, common)
    pj = np.searchsorted(idx_j, common)
    out = np.empty(common.size, np.float32)
    for lo in range(0, common.size, BLK):
        hi = min(lo + BLK, common.size)
        out[lo:hi] = (U_i[pi[lo:hi]] * U_j[pj[lo:hi]]).sum(axis=1)
    return common, out


def main() -> None:
    t0 = time.time()
    WORK.mkdir(exist_ok=True)
    ka = ka1()
    if not ka["pass_"]:
        raise SystemExit("KA-1 FAILED — abort (AGENTS.md 1375-1382)")

    # sources: 276's scan, single-line CRISPRi/KO, >= 5 held genes (PREREGISTRATION s1)
    per_src = {}
    for k, (rel, line, ds) in I.CANDIDATES.items():
        if not (ROOT / rel / "perts.txt").exists():
            continue
        perts, genes, _ = L.read_bank(rel)
        per_src[k] = dict(line=line, dataset=ds, rel=rel, perts=set(perts), genes=set(genes))
    genes_held = held_genes(per_src)
    keys = sorted(k for k, d in per_src.items()
                  if len(d["perts"] & set(genes_held)) >= 5)
    axis = L.common_axis_paths([per_src[k]["rel"] for k in keys])
    lines = [per_src[k]["line"] for k in keys]
    print(f"held genes {len(genes_held)}; sources {len(keys)}; axis {len(axis)} genes; "
          f"{len(set(lines))} lines", flush=True)

    ka_0 = ka0(axis, genes_held)
    if not ka_0["pass_"]:
        raise SystemExit("KA-0 FAILED — abort")

    idxs, Us, norms = {}, {}, {}
    for k in keys:
        t1 = time.time()
        i, u, n = load_unit(per_src[k]["rel"], axis, genes_held)
        idxs[k], Us[k], norms[k] = i, u, n
        print(f"  {k:22s} {i.size:6d} held rows  ({time.time()-t1:.1f}s)", flush=True)

    cross = [(i, j) for i in range(len(keys)) for j in range(i + 1, len(keys)) if lines[i] != lines[j]]
    within = [(i, j) for i in range(len(keys)) for j in range(i + 1, len(keys)) if lines[i] == lines[j]]
    print(f"{len(cross)} cross-line pairs, {len(within)} same-line pairs", flush=True)

    nG = len(genes_held)

    def build_C(pairs):
        C = np.full((len(pairs), nG), np.nan, np.float32)
        for p, (i, j) in enumerate(pairs):
            g, c = pair_cos(idxs[keys[i]], Us[keys[i]], idxs[keys[j]], Us[keys[j]])
            C[p, g] = c
        return C

    Ccross = build_C(cross)
    Cwithin = build_C(within)
    with np.errstate(invalid="ignore"):
        c_held = np.nanmean(Ccross, axis=0)
        n_pairs = np.isfinite(Ccross).sum(axis=0)
        c_within = np.nanmean(Cwithin, axis=0)
        n_within = np.isfinite(Cwithin).sum(axis=0)
    print(f"C_held: n {int(np.isfinite(c_held).sum())}, median {np.nanmedian(c_held):+.5f}, "
          f"mean {np.nanmean(c_held):+.5f}, share<=0 {np.nanmean(c_held[np.isfinite(c_held)] <= 0):.3f}",
          flush=True)

    # --- response magnitude: per-source scale-normalised row norm, median over sources ----------
    magM = np.full((len(keys), nG), np.nan, np.float32)
    for q, k in enumerate(keys):
        if idxs[k].size:
            magM[q, idxs[k]] = norms[k] / np.median(norms[k])
    with np.errstate(invalid="ignore"):
        mag = np.nanmedian(magM, axis=0)
        n_lines_g = np.zeros(nG, np.int32)
    seen = {}
    for q, k in enumerate(keys):
        for g in idxs[k]:
            seen.setdefault(int(g), set()).add(lines[q])
    for g, s in seen.items():
        n_lines_g[g] = len(s)

    # --- CCLE median expression -----------------------------------------------------------------
    z = np.load(ROOT / "data/external/ccle_broad_2025/ccle_cache.npz", allow_pickle=True)
    cg = {str(g): i for i, g in enumerate(z["genes"])}
    X = z["X"]
    expr = np.full(nG, np.nan)
    for i, g in enumerate(genes_held):
        j = cg.get(g)
        if j is not None:
            expr[i] = float(np.median(X[:, j]))
    del X, z
    print(f"CCLE expression available for {int(np.isfinite(expr).sum())}/{nG} held genes", flush=True)

    # --- target-shuffled null (DEFINITION-276 s4), 50 draws ------------------------------------
    rng = np.random.default_rng(SEED)
    null = np.full((N_NULL, nG), np.nan, np.float32)
    for d in range(N_NULL):
        acc = np.zeros(nG); cnt = np.zeros(nG)
        for (i, j) in cross:
            ki, kj = keys[i], keys[j]
            common = np.intersect1d(idxs[ki], idxs[kj], assume_unique=True)
            if common.size == 0:
                continue
            pi = np.searchsorted(idxs[ki], common)
            pool = np.arange(idxs[kj].size)
            pick = pool[rng.integers(0, pool.size, size=common.size)]
            clash = idxs[kj][pick] == common
            while clash.any():
                pick[clash] = pool[rng.integers(0, pool.size, size=int(clash.sum()))]
                clash = idxs[kj][pick] == common
            for lo in range(0, common.size, BLK):
                hi = min(lo + BLK, common.size)
                c = (Us[ki][pi[lo:hi]] * Us[kj][pick[lo:hi]]).sum(axis=1)
                acc[common[lo:hi]] += c; cnt[common[lo:hi]] += 1
        with np.errstate(invalid="ignore"):
            null[d] = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
        if d == 0:
            print(f"  null draw 0: median {np.nanmedian(null[0]):+.6f} ({time.time()-t0:.0f}s)", flush=True)
    with np.errstate(invalid="ignore"):
        null_mean = np.nanmean(null, axis=0)

    np.savez_compressed(
        WORK / "held279.npz",
        genes=np.asarray(genes_held, dtype=object), keys=np.asarray(keys, dtype=object),
        lines=np.asarray(lines, dtype=object), axis=np.asarray(axis, dtype=object),
        pairs_cross=np.asarray(cross), pairs_within=np.asarray(within),
        Ccross=Ccross, Cwithin=Cwithin, c_held=c_held, n_pairs=n_pairs,
        c_within=c_within, n_within=n_within, mag=mag, expr=expr,
        n_lines=n_lines_g, null_mean=null_mean,
        null_sd=np.nanstd(null, axis=0).astype(np.float32))

    summary = dict(
        n_sources=len(keys), sources=keys, lines=sorted(set(lines)), n_lines=len(set(lines)),
        n_genes_axis=len(axis), n_held=nG, n_cross_pairs=len(cross), n_within_pairs=len(within),
        pairs_per_gene=dict(min=int(n_pairs.min()), median=float(np.median(n_pairs)),
                            max=int(n_pairs.max())),
        c_held=dict(mean=float(np.nanmean(c_held)), median=float(np.nanmedian(c_held)),
                    p10=float(np.nanpercentile(c_held, 10)), p90=float(np.nanpercentile(c_held, 90)),
                    share_le0=float(np.nanmean(c_held[np.isfinite(c_held)] <= 0)),
                    n=int(np.isfinite(c_held).sum())),
        c_within=dict(median=float(np.nanmedian(c_within)), n=int(np.isfinite(c_within).sum())),
        null=dict(mean=float(np.nanmean(null_mean)), median=float(np.nanmedian(null_mean)),
                  sd_over_genes=float(np.nanstd(null_mean)), n_draws=N_NULL),
        ka1=ka, ka0=ka_0, seconds=time.time() - t0,
        peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30)
    (WORK / "build279.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k != "sources"}, indent=1))


if __name__ == "__main__":
    main()
