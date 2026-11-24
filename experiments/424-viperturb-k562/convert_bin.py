"""424 — one VIPerturb-seq genome-wide bin (Seurat v5 RDS) -> per-label summed raw counts, streamed (low memory).

`readRDS` needs 16.7 GB for binC and ~22 GB for binB, so the bins are read with `rds_stream.RDSStream`: one pass over
the gzip stream, the counts' `i` and `x` vectors spilled to disk, everything else in memory. Then, in cell chunks:

    sums[g, :]      = sum over the cells labelled g (meta.data$gene) of their raw UMI counts   (float64, exact)
    libsize[c]      = colSums(counts)[c]           (audited against meta.data$nCount_RNA: must be identical)
    nfeat[c]        = nonzeros in column c         (audited against meta.data$nFeature_RNA)
    target_count[c] = counts[own target gene, c]   (0 when absent; NaN when the target is not a feature)

Writes work/conv/<bin>/: sums.npy, groups.txt, ncells.txt, features.txt, cells.tsv.gz, info.json. The spill files are
deleted at the end. `scripts/build_donor_table.py` needs only the group sums (b = log1p(TS * P / sum P) of summed counts).

    uv run python convert_bin.py binC
"""
from __future__ import annotations

import gzip
import hashlib
import json
import resource
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from rds_stream import RDSStream, RList, S4, Spilled, Vec

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"


def names_of(x):
    n = getattr(x, "attrs", {}).get("names")
    return list(n) if n is not None else None


def named(x, key):
    return x[names_of(x).index(key)]


def dimnames(x):
    dn = x.attrs.get("dimnames")
    return (list(dn[0]) if dn is not None and dn[0] is not None else None,
            list(dn[1]) if dn is not None and len(dn) > 1 and dn[1] is not None else None)


def column(df, col):
    v = named(df, col)
    lev = getattr(v, "attrs", {}).get("levels")
    if lev is not None:                                    # factor: 1-based codes
        lev = list(lev)
        return np.asarray([lev[i - 1] if i > 0 else None for i in np.asarray(v)], dtype=object)
    return np.asarray(v, dtype=object) if isinstance(v, RList) else np.asarray(v)


def rss_gib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30   # macOS: bytes


def main(b: str, chunk_cells: int = 4000) -> None:
    t0 = time.time()
    log = lambda *a: print(f"[{b} {time.time() - t0:6.0f}s rss {rss_gib():.2f}G]", *a, flush=True)
    src = WORK / "dl" / f"genome_wide_{b}.RDS"
    out = WORK / "conv" / b
    spill = WORK / "spill" / b
    out.mkdir(parents=True, exist_ok=True)
    rd = RDSStream(src, spill, log=log)
    hdr, obj, at_eof = rd.read()
    log(f"parsed; header {hdr}; clean EOF {at_eof}; decompressed bytes {rd.bytes_read:,}")
    assert at_eof, "trailing bytes after the object"
    assert isinstance(obj, S4) and list(obj.attrs["class"]) == ["Seurat"], obj.attrs.get("class")
    assays = obj.attrs["assays"]
    rna = named(assays, "RNA")
    layers = rna.attrs["layers"]
    lnames = names_of(layers)
    log(f"assays {names_of(assays)}; RNA class {list(rna.attrs['class'])}; layers {lnames}")
    features = dimnames(rna.attrs["features"])[0]
    cells = dimnames(rna.attrs["cells"])[0]
    meta = obj.attrs["meta.data"]
    mcols = names_of(meta)
    rown = list(meta.attrs["row.names"])
    log(f"features {len(features):,}; cells {len(cells):,}; meta {len(rown):,} rows, cols {mcols}")
    assert rown == cells, "meta.data row order != RNA cells"
    assert lnames == ["counts"], f"expected one counts layer, got {lnames}"
    cnt = named(layers, "counts")
    assert list(cnt.attrs["class"]) == ["dgCMatrix"], cnt.attrs["class"]
    dim = [int(v) for v in cnt.attrs["Dim"]]
    assert dim == [len(features), len(cells)], (dim, len(features), len(cells))
    fdn, cdn = dimnames(cnt)
    assert fdn in (None, features) and cdn in (None, cells)
    p = np.asarray(cnt.attrs["p"], dtype=np.int64)
    ispill, xspill = cnt.attrs["i"], cnt.attrs["x"]
    assert isinstance(ispill, Spilled) and isinstance(xspill, Spilled)
    nnz = int(p[-1])
    assert ispill.n == nnz == xspill.n and xspill.audit["f32_exact"], (ispill.n, nnz, xspill.audit)
    log(f"counts {dim[0]:,} x {dim[1]:,}, nnz {nnz:,}; x audit {xspill.audit}")

    gene = column(meta, "gene").astype(str)
    guide = column(meta, "guide").astype(str) if "guide" in mcols else np.full(len(cells), "")
    sample = column(meta, "sample").astype(str) if "sample" in mcols else np.full(len(cells), "")
    ncount = np.asarray(column(meta, "nCount_RNA"), dtype=np.float64)
    nfeat_meta = np.asarray(column(meta, "nFeature_RNA"), dtype=np.float64)
    groups, codes, ncell = np.unique(gene, return_inverse=True, return_counts=True)
    fidx = {g: i for i, g in enumerate(features)}
    trow = np.asarray([fidx.get(g, -1) for g in gene], dtype=np.int64)
    G, F = len(groups), len(features)
    log(f"{G:,} labels; NO-TARGET cells {int(ncell[list(groups).index('NO-TARGET')]) if 'NO-TARGET' in groups else 0}; "
        f"targets not a feature: {int((trow < 0).sum())} cells")

    acc = np.zeros((G, F), dtype=np.float64)
    lib = np.zeros(len(cells), dtype=np.float64)
    nfe = np.zeros(len(cells), dtype=np.int64)
    tcount = np.zeros(len(cells), dtype=np.float64)
    with open(ispill.path, "rb") as fi, open(xspill.path, "rb") as fx:
        for c0 in range(0, len(cells), chunk_cells):
            c1 = min(c0 + chunk_cells, len(cells))
            a, z = int(p[c0]), int(p[c1])
            fi.seek(4 * a)
            fx.seek(4 * a)
            ii = np.fromfile(fi, dtype=np.int32, count=z - a)
            xx = np.fromfile(fx, dtype=np.float32, count=z - a).astype(np.float64)
            per = np.diff(p[c0:c1 + 1])
            local = np.repeat(np.arange(c1 - c0), per)
            lib[c0:c1] = np.bincount(local, weights=xx, minlength=c1 - c0)
            nfe[c0:c1] = per
            own = ii == trow[c0:c1][local]
            tcount[c0:c1] = np.bincount(local[own], weights=xx[own], minlength=c1 - c0)
            m = sp.csr_matrix((xx, (codes[c0:c1][local], ii)), shape=(G, F))
            m.sum_duplicates()
            m = m.tocoo()
            acc[m.row, m.col] += m.data
            if (c0 // chunk_cells) % 20 == 0:
                log(f"  cells {c1:,}/{len(cells):,}")
    tcount[trow < 0] = np.nan
    d_lib = float(np.abs(lib - ncount).max())
    d_nf = float(np.abs(nfe - nfeat_meta).max())
    log(f"audit: max |colSums - nCount_RNA| = {d_lib}; max |nnz/col - nFeature_RNA| = {d_nf}; "
        f"sum(acc) {acc.sum():.6e} == sum(x) {xspill.audit['sum']:.6e}")
    assert d_lib == 0.0 and d_nf == 0.0, "counts do not reproduce Seurat's own per-cell totals"
    assert abs(acc.sum() - xspill.audit["sum"]) <= 1e-9 * xspill.audit["sum"]
    np.save(out / "sums.npy", acc)
    (out / "groups.txt").write_text("\n".join(groups))
    (out / "ncells.txt").write_text("\n".join(map(str, ncell)))
    (out / "features.txt").write_text("\n".join(features))
    with gzip.open(out / "cells.tsv.gz", "wt") as fh:
        fh.write("barcode\tgene\tguide\tsample\tnCount_RNA\tlibsize\tnFeature\ttarget_count\n")
        for r in zip(cells, gene, guide, sample, ncount, lib, nfe, tcount):
            fh.write("\t".join(map(str, r)) + "\n")
    info = dict(input=src.name, input_bytes=src.stat().st_size, header=hdr, decompressed_bytes=rd.bytes_read,
                n_features=F, n_cells=len(cells), n_labels=G, nnz=nnz, x_audit=xspill.audit, meta_columns=mcols,
                layers=lnames, audit_libsize_max_abs=d_lib, audit_nfeature_max_abs=d_nf,
                features_sha256=hashlib.sha256("\n".join(features).encode()).hexdigest(),
                seconds=round(time.time() - t0), peak_rss_gib=rss_gib())
    (out / "info.json").write_text(json.dumps(info, indent=1, default=float))
    shutil.rmtree(spill)
    log(f"wrote {out}; spill removed; peak RSS {rss_gib():.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1])
