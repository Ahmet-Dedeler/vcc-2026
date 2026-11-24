"""194 — streaming accumulation of every estimator variant, in one pass per stage.

Two stages, because the knockdown weight of a cell depends on that cell's own count of
its target gene, which is only visible once the expression matrix has been read:

  stage `obs`   read the matrix, keep ONLY (per cell) its own target gene's count and its
                library size.  No accumulator, flat memory.
  stage `sums`  read the matrix again with the per-cell weight book in hand and carry one
                weighted group-sum accumulator per scheme.

`v0` reproduces the incumbent builder's accumulation (same file/block order, same
float64 sparse-indicator matmul on Orion), so its rows must come back out equal to the
shipped `data/donors/<donor>/delta.npy`.  That is the known-answer row for this whole
re-implementation: if `v0` does not reproduce, nothing else here can be read.

    uv run python accumulate.py --donor orion_hct116 --stage obs
    uv run python accumulate.py --donor orion_hct116 --stage sums --targets restricted
"""
from __future__ import annotations

import argparse
import json
import os
import re
import resource
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import h5py
import numpy as np
import pyarrow.parquet as pq
import scipy.sparse as sp

from schemes import WeightBook

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

TS = 5e4
MIN_CELLS = 20
CONTROL_PATTERNS = (r"^non.?targeting", r"^ntc", r"^control", r"^safe.?harbou?r", r"^neg")
N_EXTRA_TARGETS = 1300
# Orion carries 38,606 gene columns against gwps's 8,248, so the same accumulator count
# costs 4.7x the memory.  Each parquet file is a SINGLE row group, so pyarrow must
# decompress ~1.9 GB before `iter_batches` can slice it and the reader alone peaks at
# 5.6 GiB; the Orion row set is cut to keep the total inside the 14 GiB budget.  The
# extras are a nested prefix of the same fixed draw, so the Orion set is a subset of the
# gwps set and the two remain directly comparable.
N_EXTRA_BY_DONOR = {"k562_gwps": 1300, "orion_hct116": 500, "orion_hek293t": 500,
                    "k562_essential": 1300}
# the shared-target pool is frozen over the three PRODUCTION donors so that adding a
# further donor here cannot change a row set that has already been built
SHARED_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
SEED = 20260903
RSS_LIMIT = 13.0

DONORS = {
    "orion_hct116": dict(kind="orion", line="HCT116"),
    "orion_hek293t": dict(kind="orion", line="HEK293T"),
    "k562_gwps": dict(kind="gwps",
                      path="data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad",
                      control="control"),
    # same cell line as k562_gwps, different library: 159's same-line positive control,
    # and the denominator of its route-2 estimate of kappa^2
    "k562_essential": dict(
        kind="gwps",
        path="data/external/scperturb/ReplogleWeissman2022_K562_essential.h5ad",
        control="control"),
}


def rss_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1 << 30)


def guard() -> None:
    if rss_gib() > RSS_LIMIT:
        raise SystemExit(f"peak RSS {rss_gib():.1f} GiB over {RSS_LIMIT} GiB budget")


def is_control(name: str) -> bool:
    n = str(name).lower()
    return any(re.match(p, n) for p in CONTROL_PATTERNS)


def orion_files(line: str) -> list[Path]:
    return sorted((ROOT / "data/external/orion/data").glob(f"{line}_Batch*.parquet"),
                  key=lambda p: int(re.search(r"Batch(\d+)", p.name)[1]))


# --------------------------------------------------------------------- cell tables
def cell_table(donor: str) -> dict:
    cfg = DONORS[donor]
    if cfg["kind"] == "orion":
        z = np.load(WORK / f"prepass-{cfg['line']}.npz", allow_pickle=True)
        return dict(target=z["target"], guide=z["guide"], batch=z["batch"].astype(np.int64),
                    declared_depth=z["total_counts"].astype(np.float64),
                    pct_mt=z["pct_mt"].astype(np.float64),
                    n_genes=z["n_genes"].astype(np.float64),
                    keep=z["pass_filter"] == 1, is_ctrl=z["is_control"])
    with h5py.File(ROOT / cfg["path"], "r") as f:
        def cat(k):
            g = f["obs"][k]
            c = np.array([x.decode() if isinstance(x, bytes) else str(x)
                          for x in g["categories"][:]], dtype=object)
            return c[g["codes"][:]].astype("U")
        target = cat("perturbation")
        return dict(target=target, guide=cat("guide_id"),
                    batch=f["obs/batch"][:].astype(np.int64),
                    declared_depth=f["obs/ncounts"][:].astype(np.float64),
                    pct_mt=(100.0 * f["obs/percent_mito"][:]).astype(np.float64),
                    n_genes=f["obs/ngenes"][:].astype(np.float64),
                    keep=np.ones(len(target), bool),
                    is_ctrl=target == cfg["control"])


def donor_axis(donor: str) -> tuple[np.ndarray, np.ndarray]:
    """Gene symbols and the incumbent control expression fraction (read-only)."""
    d = ROOT / "data/donors" / donor
    genes = np.array((d / "genes.txt").read_text().split("\n"), dtype=object)
    b = np.load(d / "b_ctrl.npy").astype(np.float64)
    assert len(genes) == len(b), (len(genes), len(b))
    return genes, np.expm1(b) / TS


def target_column(genes: np.ndarray, p_ctrl: np.ndarray, targets: np.ndarray) -> np.ndarray:
    """Per-cell column of the cell's own target gene; -1 where the gene is not measured."""
    tok: dict[str, int] = {}
    for i, g in enumerate(genes):
        g = str(g)
        if not g:
            continue
        j = tok.get(g)
        if j is None or p_ctrl[i] > p_ctrl[j]:
            tok[g] = i
    uniq, inv = np.unique(targets, return_inverse=True)
    lut = np.array([tok.get(str(u), -1) for u in uniq], dtype=np.int64)
    return lut[inv]


# --------------------------------------------------------------------- stage: obs
def stage_obs(donor: str) -> None:
    cfg = DONORS[donor]
    cells = cell_table(donor)
    genes, p_ctrl = donor_axis(donor)
    n_genes, n = len(genes), len(cells["target"])
    col = target_column(genes, p_ctrl, cells["target"])
    obs = np.full(n, np.nan, dtype=np.float32)
    depth = np.zeros(n, dtype=np.float64)

    t0 = time.time()
    if cfg["kind"] == "orion":
        files = orion_files(cfg["line"])
        off = 0
        for i, f in enumerate(files, 1):
            for lo, X in orion_blocks(f, n_genes):
                m = X.shape[0]
                s = off + lo
                depth[s:s + m] = np.asarray(X.sum(1)).ravel()
                obs[s:s + m] = _gather(X, col[s:s + m])
                del X
            off += _orion_rows(f)
            if i % 20 == 0 or i == len(files):
                print(f"  {i}/{len(files)}  {time.time()-t0:.0f}s  rss {rss_gib():.1f}G",
                      flush=True)
                guard()
        assert off == n, (off, n)
    else:
        with h5py.File(ROOT / cfg["path"], "r") as f:
            X = f["X"]
            for lo in range(0, n, 16_384):
                hi = min(lo + 16_384, n)
                blk = np.asarray(X[lo:hi], dtype=np.float64)
                depth[lo:hi] = blk.sum(1)
                c = col[lo:hi]
                have = np.flatnonzero(c >= 0)
                obs[lo + have] = blk[have, c[have]]
                if (lo // 16_384) % 30 == 0:
                    print(f"  {hi:,}/{n:,}  {time.time()-t0:.0f}s  rss {rss_gib():.1f}G",
                          flush=True)
                    guard()

    exp = np.where(col >= 0, p_ctrl[np.maximum(col, 0)] * depth, np.nan)
    live = np.isfinite(exp) & ~cells["is_ctrl"] & cells["keep"]
    summary = dict(
        donor=donor, n_cells=int(n),
        targets_with_own_gene=int(len(np.unique(col[col >= 0]))),
        n_targets=int(len(np.unique(cells["target"][~cells["is_ctrl"]]))),
        cells_with_readout=float(live.mean()),
        matrix_vs_declared_depth=float(np.median(
            depth[cells["keep"]] / np.maximum(cells["declared_depth"][cells["keep"]], 1))),
        expected_counts=dict(p25=float(np.percentile(exp[live], 25)),
                             median=float(np.median(exp[live])),
                             p75=float(np.percentile(exp[live], 75)),
                             frac_ge_0p5=float((exp[live] >= 0.5).mean())),
        elapsed_s=round(time.time() - t0, 1), peak_rss_gib=round(rss_gib(), 2))
    np.savez_compressed(WORK / f"obs-{donor}.npz", obs=obs, exp=exp.astype(np.float32),
                        col=col, depth=depth.astype(np.float32),
                        summary=json.dumps(summary))
    (WORK / f"obs-{donor}.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


ORION_CHUNK = 4_096


def orion_blocks(path: Path, n_genes: int):
    """Yield (offset, csr_matrix) over one Orion parquet file, in file row order.

    Reads through Arrow record batches and takes the list columns' flat value/offset
    buffers directly.  Converting the list columns to pandas first (what the incumbent
    builder does) materialises one Python-level numpy array per cell and pushed peak RSS
    to 27 GiB on the first HCT116 run of this experiment; this path holds ~3 GiB.
    """
    pf = pq.ParquetFile(path)
    off = 0
    for rb in pf.iter_batches(batch_size=ORION_CHUNK,
                              columns=["gene_token_id", "gene_expression"]):
        a = rb.column("gene_token_id")
        b = rb.column("gene_expression")
        offsets = np.asarray(a.offsets, dtype=np.int64)
        indptr = offsets - offsets[0]
        idx = a.flatten().to_numpy(zero_copy_only=False).astype(np.int32, copy=False)
        val = b.flatten().to_numpy(zero_copy_only=False).astype(np.float64, copy=False)
        m = rb.num_rows
        assert len(idx) == indptr[-1] == len(val), (len(idx), indptr[-1], len(val))
        yield off, sp.csr_matrix((val, idx, indptr), shape=(m, n_genes))
        off += m


def _orion_rows(path: Path) -> int:
    return pq.ParquetFile(path).metadata.num_rows


def _gather(X: sp.csr_matrix, col: np.ndarray, chunk: int = 2_048) -> np.ndarray:
    """X[i, col[i]] for every row; NaN where col < 0.

    After `sort_indices` the flattened key `row * n_genes + index` is ascending, so one
    `searchsorted` answers a whole row block.  Done in row blocks so the key array stays
    ~10^7 elements instead of ~10^8.
    """
    X.sort_indices()
    m, K = X.shape
    out = np.full(m, np.nan)
    for lo in range(0, m, chunk):
        hi = min(lo + chunk, m)
        a, b = int(X.indptr[lo]), int(X.indptr[hi])
        if b <= a:
            continue
        loc = np.flatnonzero(col[lo:hi] >= 0)
        if not len(loc):
            continue
        rowid = np.repeat(np.arange(hi - lo, dtype=np.int64),
                          np.diff(X.indptr[lo:hi + 1]))
        key = rowid * K + X.indices[a:b].astype(np.int64)
        want = loc.astype(np.int64) * K + col[lo:hi][loc]
        pos = np.searchsorted(key, want)
        ok = (pos < len(key)) & (key[np.minimum(pos, len(key) - 1)] == want)
        v = np.zeros(len(loc))
        v[ok] = X.data[a:b][pos[ok]]
        out[lo:hi][loc] = v
    return out


# -------------------------------------------------------------------- stage: sums
def choose_targets(donor: str, mode: str) -> set[str] | None:
    if mode == "all":
        return None
    panel = [l.split(",")[0] for l in
             (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
    shared: set[str] | None = None
    for d in SHARED_DONORS:
        s = set((ROOT / "data/donors" / d / "perts.txt").read_text().split("\n"))
        shared = s if shared is None else (shared & s)
    mine = set((ROOT / "data/donors" / donor / "perts.txt").read_text().split("\n"))
    keep = {g for g in panel if g in mine}
    # the extras pool must NOT depend on the donor, or two donors draw disjoint extras
    # and their tables end up sharing almost no rows.  `shared - keep` was donor-dependent
    # only through `keep`, and for the three production donors `shared & panel` is already
    # inside `keep`, so this is a no-op for them and a fix for any donor that covers no
    # panel target (k562_essential covers 0/300).
    pool = sorted(shared - set(panel))
    rng = np.random.default_rng(SEED)
    extra = rng.choice(np.array(pool), size=min(N_EXTRA_TARGETS, len(pool)),
                       replace=False)[:N_EXTRA_BY_DONOR[donor]]
    return keep | {str(x) for x in extra}


def stage_sums(donor: str, mode: str) -> None:
    cfg = DONORS[donor]
    cells = cell_table(donor)
    genes, p_ctrl = donor_axis(donor)
    n_genes, n = len(genes), len(cells["target"])
    z = np.load(WORK / f"obs-{donor}.npz")
    obs, exp, depth = (z["obs"].astype(np.float64), z["exp"].astype(np.float64),
                       z["depth"].astype(np.float64))

    keep = cells["keep"]
    ctrl = cells["is_ctrl"] & keep
    pert = (~cells["is_ctrl"]) & keep

    tnames, tcounts = np.unique(cells["target"][pert], return_counts=True)
    eligible = {str(a) for a, b in zip(tnames, tcounts) if b >= MIN_CELLS}
    wanted = choose_targets(donor, mode)
    if wanted is not None:
        eligible &= wanted
    names = np.array(sorted(eligible))
    ctrl_code = len(names)
    tcode = np.full(n, -1, dtype=np.int64)
    if len(names):
        pos = np.searchsorted(names, cells["target"])
        ok = pert & (pos < len(names)) & (names[np.minimum(pos, len(names) - 1)] ==
                                          cells["target"])
        tcode[ok] = pos[ok]
    tcode[ctrl] = ctrl_code
    n_groups = len(names) + 1
    print(f"{donor}/{mode}: {len(names)} targets, {int(ctrl.sum()):,} control cells, "
          f"{int((tcode >= 0).sum()):,} accumulated cells", flush=True)

    book = WeightBook(depth=depth, pct_mt=cells["pct_mt"], n_genes=cells["n_genes"],
                      batch=cells["batch"], tcode=tcode, ctrl_code=ctrl_code,
                      obs_target=obs, exp_target=exp, seed=SEED)
    schemes = book.names()
    W = np.stack([book.w[s] for s in schemes])
    assert np.allclose(W[schemes.index("v0")][tcode >= 0], 1.0)
    print(f"schemes: {schemes}", flush=True)
    audit = book.audit
    print(json.dumps(audit, indent=2), flush=True)

    # NEGR1 lands in the incumbent control pool via the `^neg` regex; carry a clean
    # control vector alongside so the size of that defect can be priced.
    ctrl_clean = ctrl & (cells["target"] != "NEGR1")
    negr1_cells = int((ctrl & (cells["target"] == "NEGR1")).sum())

    # everything that reads the per-cell tables is done here; free them before the
    # accumulators are allocated, so peak RSS is accumulators + one buffered block
    live = tcode >= 0
    neff = np.zeros((len(schemes), n_groups))
    for k in range(len(schemes)):
        a = np.bincount(tcode[live], weights=W[k, live] * depth[live], minlength=n_groups)
        b = np.bincount(tcode[live], weights=W[k, live] ** 2 * depth[live],
                        minlength=n_groups)
        neff[k] = np.where(b > 0, a * a / np.maximum(b, 1e-300), 0.0)
    cellcount = np.bincount(tcode[live], minlength=n_groups)
    have = live & (tcode != ctrl_code) & np.isfinite(exp)
    kd_target = np.stack([
        np.where(np.bincount(tcode[have], weights=exp[have], minlength=n_groups) > 0,
                 np.bincount(tcode[have], weights=obs[have], minlength=n_groups) /
                 np.maximum(np.bincount(tcode[have], weights=exp[have],
                                        minlength=n_groups), 1e-12), np.nan),
        np.bincount(tcode[have], weights=exp[have], minlength=n_groups),
        np.bincount(tcode[have], minlength=n_groups).astype(float),
        getattr(book, "_kd_median_expected", np.zeros(n_groups))])
    book.w.clear()
    for k in ("guide", "target", "declared_depth", "pct_mt", "n_genes", "batch",
              "is_ctrl", "keep"):
        cells.pop(k, None)
    del obs, exp, depth, z, keep, pert, live, have, book
    import gc
    gc.collect()

    acc = np.zeros((len(schemes), n_groups, n_genes), dtype=np.float64)
    wsum = np.zeros((len(schemes), n_groups), dtype=np.float64)
    acc_ctrl_clean = np.zeros(n_genes, dtype=np.float64)
    print(f"accumulators {acc.nbytes / 2**30:.2f} GiB, rss {rss_gib():.1f}G", flush=True)

    t0 = time.time()
    if cfg["kind"] == "orion":
        files = orion_files(cfg["line"])
        off = 0
        # Only ~13% of Orion cells belong to a restricted target or to the control, and
        # the group-sum product is dense in the gene axis, so the dense conversion --
        # not the read -- dominates.  Buffering the KEPT cells across parquet chunks and
        # flushing at ~64k cuts that conversion cost by an order of magnitude.
        buf: list = []
        buf_n = 0

        def flush():
            nonlocal buf, buf_n
            if not buf:
                return
            X = sp.vstack([b[0] for b in buf], format="csr")
            codes = np.concatenate([b[1] for b in buf])
            Wb = np.concatenate([b[2] for b in buf], axis=1)
            cln = np.concatenate([b[3] for b in buf])
            _accumulate(acc, wsum, X, codes, Wb, n_groups)
            if cln.any():
                acc_ctrl_clean[:] += np.asarray(X[cln].sum(0)).ravel()
            buf, buf_n = [], 0

        for i, f in enumerate(files, 1):
            for lo, X in orion_blocks(f, n_genes):
                m = X.shape[0]
                s = off + lo
                sel = np.flatnonzero(tcode[s:s + m] >= 0)
                if len(sel):
                    buf.append((X[sel], tcode[s:s + m][sel], W[:, s:s + m][:, sel],
                                ctrl_clean[s:s + m][sel]))
                    buf_n += len(sel)
                del X
                if buf_n >= 6_000:
                    flush()
            off += _orion_rows(f)
            if i % 10 == 0 or i == len(files):
                print(f"  {i}/{len(files)}  {time.time()-t0:.0f}s  rss {rss_gib():.1f}G",
                      flush=True)
                guard()
        flush()
        assert off == n
    else:
        with h5py.File(ROOT / cfg["path"], "r") as f:
            Xh = f["X"]
            for lo in range(0, n, 16_384):
                hi = min(lo + 16_384, n)
                blk = np.asarray(Xh[lo:hi], dtype=np.float64)
                _accumulate(acc, wsum, blk, tcode[lo:hi], W[:, lo:hi], n_groups)
                cc = ctrl_clean[lo:hi]
                if cc.any():
                    acc_ctrl_clean += blk[cc].sum(0)
                if (lo // 16_384) % 30 == 0:
                    print(f"  {hi:,}/{n:,}  {time.time()-t0:.0f}s  rss {rss_gib():.1f}G",
                          flush=True)
                    guard()

    out = WORK / f"sums-{donor}-{mode}"
    out.mkdir(parents=True, exist_ok=True)
    for k, s in enumerate(schemes):
        np.save(out / f"{s}.npy", acc[k])
    np.save(out / "_ctrl_clean.npy", acc_ctrl_clean)
    np.save(out / "_wsum.npy", wsum)
    np.save(out / "_cells.npy", cellcount)
    # effective sample size per (scheme, group): (sum w*d)^2 / sum w^2*d, used by the
    # Poisson log1p debias correction
    np.save(out / "_neff.npy", neff)
    # per-target knockdown diagnostic: pooled residual expression of the target's own
    # gene, sum(observed) / sum(expected-under-control) -- the target-level version of
    # the readout the cell weights use.
    np.save(out / "_kd_target.npy", kd_target)
    (out / "meta.json").write_text(json.dumps(dict(
        donor=donor, mode=mode, schemes=schemes, names=names.tolist(),
        ctrl_code=ctrl_code, n_genes=int(n_genes), n_control_cells=int(ctrl.sum()),
        n_control_cells_clean=int(ctrl_clean.sum()), negr1_in_control=negr1_cells,
        audit=audit,
        elapsed_s=round(time.time() - t0, 1), peak_rss_gib=round(rss_gib(), 2)), indent=2))
    print(f"wrote {out}  elapsed {time.time()-t0:.0f}s  peak rss {rss_gib():.2f} GiB")


def _accumulate(acc, wsum, X, codes, Wblk, n_groups) -> None:
    """acc[k, g] += sum_{c in g} W[k,c] * x_c, for every scheme k."""
    live = codes >= 0
    if not live.any():
        return
    if live.all():
        rows, Xp = np.arange(len(codes)), X
    else:
        rows, Xp = np.flatnonzero(live), X[np.flatnonzero(live)]
    g = codes[rows]
    sparse = sp.issparse(Xp)
    for k in range(acc.shape[0]):
        w = Wblk[k, rows]
        nz = np.flatnonzero(w != 0.0)
        if not len(nz):
            continue
        ug, inv = np.unique(g[nz], return_inverse=True)
        ind = sp.csr_matrix((w[nz], (inv, nz)), shape=(len(ug), len(rows)))
        prod = (ind @ Xp).toarray() if sparse else ind @ Xp
        # row-by-row so no (len(ug), n_genes) fancy-index temporary is materialised
        for j, gi in enumerate(ug):
            acc[k][gi] += prod[j]
        del prod, ind
        wsum[k] += np.bincount(g[nz], weights=w[nz], minlength=n_groups)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", required=True, choices=list(DONORS))
    ap.add_argument("--stage", required=True, choices=["obs", "sums"])
    ap.add_argument("--targets", default="restricted", choices=["restricted", "all"])
    args = ap.parse_args()
    (stage_obs(args.donor) if args.stage == "obs"
     else stage_sums(args.donor, args.targets))


if __name__ == "__main__":
    main()
