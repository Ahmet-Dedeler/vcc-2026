"""161 SS4 — per-gene CONTROL-CELL statistics per cell line, from cell-level data.

WHY THIS FILE EXISTS
--------------------
A donor table stores only `b_ctrl`, the control *pseudobulk*.  159's claim is that the
information production is missing is CONTEXT, and the two context scalars a donor table
cannot express are

    detection rate   fraction of control cells in which the gene is observed at all
    dispersion       how variable the gene is ACROSS control cells, beyond its mean

Both are per-(line, gene) scalars, both are computable in the challenge from the 18,400
control cells we are handed, and neither is a function of the pseudobulk.  160 showed the
control cells' dominant *direction* is unusable; this asks a different question — whether
their per-gene second moment says which genes the line's response will load on.

WHAT IS COMPUTED, PER LINE, ON THAT LINE'S OWN DONOR GENE AXIS
--------------------------------------------------------------
    det       mean over control cells of 1[count > 0]
    mlog      mean over control cells of log1p(TS * count / libsize)      (scorer's space)
    vlog      variance of the same, ddof=1
    fano      var(cpm) / mean(cpm)   -- overdispersion, scale-free in the count space
    b_ctrl_c  log1p(TS * sum(count) / sum(libsize)), the pseudobulk rebuilt HERE

`b_ctrl_c` is the known-answer row of this script: it is computed from the same cells the
donor builder used and must reproduce the shipped `b_ctrl.npy` closely.  If it does not, the
control cell set identified here is not the one the donor table was built from and every
statistic below is about the wrong cells.

Streaming: dense h5ad is read in row slabs and only control rows are kept; Orion is
accumulated straight off the parquet sparse lists, exactly as `build_orion_donor.py` does.
Peak RSS is one slab, ~1 GB.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import time

import numpy as np

from common import REPO, TS, WORK, load_table

SLAB = 20_000
CONTROL_PATTERNS = (r"^non.?targeting$", r"^ntc$", r"^control$", r"^safe.?harbou?r$", r"^neg")

# line -> (kind, path, obs column holding the perturbation label)
SOURCES = {
    # K562 control cells come from the ESSENTIAL library, not gwps: gwps is a
    # (1.99M x 8248) dense array chunked (3886, 33), so pulling its scattered control rows
    # costs a full 65 GB of chunk reads.  Same cell line, same lab, same chemistry; the
    # statistic wanted here is a property of K562's basal state, not of the library.
    "K562":    ("h5ad", "data/external/scperturb/ReplogleWeissman2022_K562_essential.h5ad",
                "perturbation"),
    "RPE1":    ("h5ad", "data/external/scperturb/ReplogleWeissman2022_rpe1.h5ad",
                "perturbation"),
    "JURKAT":  ("h5ad", "data/external/nadig2025/jurkat_raw_singlecell.h5ad", "gene"),
    "HCT116":  ("orion", "HCT116", None),
    "HEK293T": ("orion", "HEK293T", None),
    "H1":      ("h5ad_harness", "data/harness/h1-test.h5ad", "target_gene"),
}

# which donor table defines each line's gene axis
AXIS = {"K562": "k562_gwps", "RPE1": "rpe1", "JURKAT": "nadig_jurkat",
        "HCT116": "orion_hct116", "HEK293T": "orion_hek293t", "H1": "h1_2025"}


def is_control(name: str) -> bool:
    n = str(name).lower()
    return any(re.match(p, n) for p in CONTROL_PATTERNS)


class Acc:
    """Streaming per-gene accumulator.  Everything is one pass, O(n_genes) memory."""

    def __init__(self, n_genes: int):
        self.n = 0
        self.nz = np.zeros(n_genes)          # cells with a nonzero count
        self.s_cnt = np.zeros(n_genes)       # sum of raw counts
        self.s_cpm = np.zeros(n_genes)       # sum of per-cell CPM-like (TS-scaled) values
        self.s_cpm2 = np.zeros(n_genes)
        self.s_log = np.zeros(n_genes)       # sum of log1p(TS * p)
        self.s_log2 = np.zeros(n_genes)
        self.s_lib = 0.0

    def add_dense(self, block: np.ndarray) -> None:
        lib = block.sum(1, keepdims=True)
        ok = lib.ravel() > 0
        if not ok.all():
            block, lib = block[ok], lib[ok]
        if not len(block):
            return
        q = (TS * block) / lib                       # TS-scaled relative abundance
        lg = np.log1p(q)
        self.n += len(block)
        self.nz += (block > 0).sum(0)
        self.s_cnt += block.sum(0)
        self.s_cpm += q.sum(0)
        self.s_cpm2 += (q * q).sum(0)
        self.s_log += lg.sum(0)
        self.s_log2 += (lg * lg).sum(0)
        self.s_lib += float(lib.sum())

    def add_csr(self, X) -> None:
        """Same accumulation without densifying: rows are short sparse lists."""
        import scipy.sparse as sp
        lib = np.asarray(X.sum(1)).ravel()
        keep = lib > 0
        if not keep.all():
            X, lib = X[keep], lib[keep]
        if X.shape[0] == 0:
            return
        inv = sp.diags(TS / lib)
        Q = (inv @ X).tocsr()                        # TS * count / libsize, still sparse
        L = Q.copy()
        L.data = np.log1p(L.data)                    # log1p(0) = 0, so sparsity is preserved
        self.n += X.shape[0]
        self.nz += np.asarray((X > 0).sum(0)).ravel()
        self.s_cnt += np.asarray(X.sum(0)).ravel()
        self.s_cpm += np.asarray(Q.sum(0)).ravel()
        Q2 = Q.copy(); Q2.data = Q2.data ** 2
        self.s_cpm2 += np.asarray(Q2.sum(0)).ravel()
        self.s_log += np.asarray(L.sum(0)).ravel()
        L2 = L; L2.data = L2.data ** 2
        self.s_log2 += np.asarray(L2.sum(0)).ravel()
        self.s_lib += float(lib.sum())

    def finish(self) -> dict:
        n = max(self.n, 1)
        mlog = self.s_log / n
        vlog = np.maximum(self.s_log2 / n - mlog ** 2, 0.0) * (n / max(n - 1, 1))
        mcpm = self.s_cpm / n
        vcpm = np.maximum(self.s_cpm2 / n - mcpm ** 2, 0.0) * (n / max(n - 1, 1))
        return dict(n_cells=self.n, det=self.nz / n, mlog=mlog, vlog=vlog,
                    fano=vcpm / np.maximum(mcpm, 1e-12),
                    b_ctrl_c=np.log1p(TS * self.s_cnt / max(self.s_lib, 1.0)))


def from_h5ad(path: str, obs_col: str, gene_names: list[str], harness: bool) -> dict:
    import h5py
    p = os.path.join(REPO, path)
    with h5py.File(p, "r") as f:
        o = f["obs"][obs_col]
        if isinstance(o, h5py.Group):                       # anndata categorical
            cats = [c.decode() if isinstance(c, bytes) else str(c)
                    for c in o["categories"][:]]
            codes = o["codes"][:]
        else:                                               # legacy 0.1.0 __categories
            codes = o[:]
            cats = [c.decode() if isinstance(c, bytes) else str(c)
                    for c in f["obs"]["__categories"][obs_col][:]]
        ctrl_cat = np.asarray([is_control(c) for c in cats])
        mask = ctrl_cat[codes]
        print(f"   control labels: {[c for c, k in zip(cats, ctrl_cat) if k]}"
              f"  -> {int(mask.sum()):,} of {len(mask):,} cells")
        assert mask.sum() >= 500, "too few control cells identified"

        # gene axis of the file -> position in the donor table's axis.  Three storage
        # conventions occur across these files and the wrong one silently maps 0 genes:
        # plain string array, anndata-0.2 categorical group, and legacy-0.1 __categories.
        def var_col(name):
            v = f["var"][name]
            if isinstance(v, h5py.Group):                    # 0.2 categorical
                cats = [c.decode() if isinstance(c, bytes) else str(c)
                        for c in v["categories"][:]]
                return [cats[i] for i in v["codes"][:]]
            arr = v[:]
            if arr.dtype.kind in "iu" and "__categories" in f["var"]:   # 0.1 categorical
                cats = [c.decode() if isinstance(c, bytes) else str(c)
                        for c in f["var"]["__categories"][name][:]]
                return [cats[i] for i in arr]
            return [g.decode() if isinstance(g, bytes) else str(g) for g in arr]

        vi = f["var"].attrs.get("_index", "_index")
        fg = var_col("gene_name") if "gene_name" in f["var"] else var_col(vi)
        pos = {g: i for i, g in enumerate(fg)}
        col = np.asarray([pos.get(g, -1) for g in gene_names])
        hit = col >= 0
        print(f"   gene axis: {hit.sum():,}/{len(gene_names):,} of the donor axis present")

        X = f["X"]
        acc = Acc(len(fg))
        if isinstance(X, h5py.Group):                       # sparse CSR on disk
            import scipy.sparse as sp
            indptr = X["indptr"][:]
            nrow = len(indptr) - 1
            for s in range(0, nrow, SLAB):
                e = min(s + SLAB, nrow)
                sel = np.where(mask[s:e])[0]
                if not len(sel):
                    continue
                lo, hi = indptr[s], indptr[e]
                sub = sp.csr_matrix((X["data"][lo:hi], X["indices"][lo:hi],
                                     indptr[s:e + 1] - lo), shape=(e - s, len(fg)))
                acc.add_csr(sub[sel])
        else:
            nrow = X.shape[0]
            for s in range(0, nrow, SLAB):
                e = min(s + SLAB, nrow)
                sel = np.where(mask[s:e])[0]
                if not len(sel):
                    continue
                acc.add_dense(np.asarray(X[s:e][sel], dtype=np.float64))
    out = acc.finish()
    return {k: (v if np.isscalar(v) else _project(v, col, hit)) for k, v in out.items()}


def _project(v: np.ndarray, col: np.ndarray, hit: np.ndarray) -> np.ndarray:
    o = np.full(len(col), np.nan)
    o[hit] = v[col[hit]]
    return o


def from_orion(line: str, gene_names: list[str]) -> dict:
    import pandas as pd
    import pyarrow.parquet as pq
    import scipy.sparse as sp
    from pathlib import Path

    gm = pd.read_parquet(os.path.join(REPO, "data/external/orion/metadata/gene_metadata.parquet"))
    gm = gm.sort_values("gene_token_id")
    n_tok = int(gm["gene_token_id"].max()) + 1
    names = np.array([""] * n_tok, dtype=object)
    names[gm["gene_token_id"].to_numpy()] = gm["gene_name"].to_numpy()
    pos = {g: i for i, g in enumerate(names)}
    col = np.asarray([pos.get(g, -1) for g in gene_names])
    hit = col >= 0
    print(f"   gene axis: {hit.sum():,}/{len(gene_names):,} of the donor axis present")

    files = sorted(Path(os.path.join(REPO, "data/external/orion/data")).glob(f"{line}_Batch*.parquet"),
                   key=lambda p: int(re.search(r"Batch(\d+)", p.name)[1]))
    assert files, f"no parquet for {line}"
    acc = Acc(n_tok)
    for i, fp in enumerate(files, 1):
        t = pq.read_table(fp, columns=["gene_token_id", "gene_expression", "gene_target",
                                       "pass_guide_filter"]).to_pandas()
        t = t[t["pass_guide_filter"] == 1]
        if t.empty:
            continue
        m = np.fromiter((is_control(x) for x in t["gene_target"].to_numpy()), bool, len(t))
        if not m.any():
            continue
        t = t[m]
        lengths = t["gene_token_id"].map(len).to_numpy()
        indptr = np.zeros(len(t) + 1, dtype=np.int64)
        np.cumsum(lengths, out=indptr[1:])
        idx = np.concatenate(t["gene_token_id"].to_numpy()).astype(np.int32)
        val = np.concatenate(t["gene_expression"].to_numpy()).astype(np.float64)
        acc.add_csr(sp.csr_matrix((val, idx, indptr), shape=(len(t), n_tok)))
        if i % 25 == 0:
            print(f"   batch {i}/{len(files)}  {acc.n:,} control cells", flush=True)
    print(f"   {acc.n:,} control cells over {len(files)} batches")
    out = acc.finish()
    return {k: (v if np.isscalar(v) else _project(v, col, hit)) for k, v in out.items()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("lines", nargs="*", default=list(SOURCES))
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    summary = {}
    for line in (a.lines or list(SOURCES)):
        p = os.path.join(WORK, f"ctrl-{line}.npz")
        if os.path.exists(p):
            print(f"[{line}] cached, skipping")
            continue
        kind, src, obs_col = SOURCES[line]
        tab = load_table(AXIS[line])
        genes = tab["genes"]
        t0 = time.time()
        print(f"[{line}] {kind} {src}  ({len(genes):,} genes on the {AXIS[line]} axis)")
        if kind == "orion":
            st = from_orion(src, genes)
        else:
            st = from_h5ad(src, obs_col, genes, kind == "h5ad_harness")

        # ---- known-answer row: rebuild the pseudobulk and compare to the shipped b_ctrl
        ok = np.isfinite(st["b_ctrl_c"])
        r = float(np.corrcoef(st["b_ctrl_c"][ok], tab["b"][ok])[0, 1])
        mad = float(np.abs(st["b_ctrl_c"][ok] - tab["b"][ok]).mean())
        print(f"   KNOWN-ANSWER b_ctrl rebuilt from these cells vs shipped b_ctrl.npy: "
              f"pearson {r:.6f}, mean|diff| {mad:.4f}  ({int(ok.sum()):,} genes)")
        assert r > 0.97, (f"{line}: rebuilt control pseudobulk does not match the donor "
                          f"table's (r={r:.4f}) -- wrong control cells")
        np.savez_compressed(p, genes=np.asarray(genes),
                            **{k: np.asarray(v) for k, v in st.items()})
        summary[line] = dict(n_cells=int(st["n_cells"]), pearson_b=r, mad_b=mad,
                             n_genes=int(ok.sum()), secs=round(time.time() - t0, 1))
        print(f"   wrote {p}  in {time.time() - t0:.0f}s\n", flush=True)
    if summary:
        with open(os.path.join(WORK, "ctrlstats.json"), "w") as fh:
            json.dump(summary, fh, indent=1)
        print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
