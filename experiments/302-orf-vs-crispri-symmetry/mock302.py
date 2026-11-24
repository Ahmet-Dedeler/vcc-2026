"""302 — a synthetic atlas with the reserve's exact schema, used to exercise `access302.py` and
`analyse302.py` end to end BEFORE the single access.  Reads only `var/gene_symbol` and
`obs/perturbation/categories` from the reserve (metadata, free); no `X`, no `obs` codes.

    PYTHONDONTWRITEBYTECODE=1 uv run python mock302.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib302 as L  # noqa: E402

N_SYM = 40
N_PER = 130
N_CTRL = 400
NNZ = 300
SPIKE = 60


def main() -> None:
    with h5py.File(L.ATLAS, "r") as f:
        var_sym = np.asarray(f["var/gene_symbol"][:]).astype(str)
        cats = np.asarray(f["obs/perturbation/categories"][:]).astype(str)
    nG = len(var_sym)
    vpos = {g: i for i, g in enumerate(var_sym)}
    code2sym = L.read_orf_map()
    perts, h1_genes, mask0 = L.read_h1_bank()
    h1_set = set(perts)
    axis_set = {g for g, m in zip(h1_genes, mask0) if m and g in vpos}

    sym_codes: dict[str, list[int]] = {}
    for c, s in code2sym.items():
        if s in L.CTRL_SYMBOLS or s not in h1_set:
            continue
        sym_codes.setdefault(s, []).append(c)
    pick = [s for s in sorted(sym_codes) if s in axis_set and len(sym_codes[s]) >= 2][:N_SYM]
    ctrl = {s: c for c, s in code2sym.items() if s in L.CTRL_SYMBOLS}
    print(f"mock: {len(pick)} symbols, {N_PER} cells each (2 isoform codes), "
          f"{N_CTRL} GFP + {N_CTRL} mCherry")

    rng = np.random.default_rng(7)
    codes, own_gene = [], []
    for s in pick:
        cs = sym_codes[s][:2]
        for j in range(N_PER):
            codes.append(cs[j % 2]); own_gene.append(vpos[s])
    for _ in range(N_CTRL):
        codes.append(ctrl["GFP"]); own_gene.append(-1)
    for _ in range(N_CTRL):
        codes.append(ctrl["mCherry"]); own_gene.append(-1)
    codes = np.asarray(codes, np.int16)
    own_gene = np.asarray(own_gene)
    n = len(codes)

    # a shared "expression programme" so residuals are not pure noise
    base = rng.random(nG) < 0.02
    pool = np.flatnonzero(base)
    indptr = np.zeros(n + 1, np.int64)
    ind_all, dat_all = [], []
    for i in range(n):
        cols = rng.choice(pool, size=min(NNZ, pool.size), replace=False)
        vals = rng.poisson(3.0, size=cols.size).astype(np.float32) + 1
        if own_gene[i] >= 0:
            cols = np.append(cols, own_gene[i])
            vals = np.append(vals, np.float32(SPIKE))
        o = np.argsort(cols)
        cols, vals = cols[o], vals[o]
        u, first = np.unique(cols, return_index=True)
        ind_all.append(u.astype(np.int64)); dat_all.append(vals[first])
        indptr[i + 1] = indptr[i] + u.size
    indices = np.concatenate(ind_all)
    data = np.concatenate(dat_all)
    ncounts = np.add.reduceat(data.astype(np.int64), indptr[:-1]).astype(np.int64)

    out = L.WORK / "mock_atlas.h5ad"
    with h5py.File(out, "w") as g:
        X = g.create_group("X")
        X.attrs["encoding-type"] = "csr_matrix"
        X.attrs["encoding-version"] = "0.1.0"
        X.attrs["shape"] = np.asarray([n, nG], np.int64)
        X.create_dataset("data", data=data)
        X.create_dataset("indices", data=indices)
        X.create_dataset("indptr", data=indptr)
        o = g.create_group("obs")
        o.attrs["_index"] = "cell_barcode"
        o.attrs["encoding-type"] = "dataframe"
        o.attrs["encoding-version"] = "0.2.0"
        o.attrs["column-order"] = np.asarray(["perturbation", "ncounts"], dtype=h5py.string_dtype())
        pg = o.create_group("perturbation")
        pg.attrs["encoding-type"] = "categorical"
        pg.create_dataset("categories", data=np.asarray(cats, dtype=object),
                          dtype=h5py.string_dtype())
        pg.create_dataset("codes", data=codes)
        o.create_dataset("ncounts", data=ncounts)
        o.create_dataset("cell_barcode", data=np.asarray([f"c{i}" for i in range(n)], dtype=object),
                         dtype=h5py.string_dtype())
        vgrp = g.create_group("var")
        vgrp.attrs["_index"] = "gene_symbol"
        vgrp.attrs["encoding-type"] = "dataframe"
        vgrp.attrs["encoding-version"] = "0.2.0"
        vgrp.attrs["column-order"] = np.asarray([], dtype=h5py.string_dtype())
        vgrp.create_dataset("gene_symbol", data=np.asarray(var_sym, dtype=object),
                            dtype=h5py.string_dtype())
    print(f"wrote {out} ({out.stat().st_size / 2**20:.1f} MiB), {n} cells, nnz {len(data)}")


if __name__ == "__main__":
    sys.exit(main())
