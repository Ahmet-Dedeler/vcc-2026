"""595 — build two mock D/E/F bundles from REAL control cells of non-A/B/C lines, in the exact bundle format.

Every rehearsal so far used A/B/C's own controls (relabelled clones, R6/R7a included). These bundles hold other lines'
control cells, on the official 18,533-gene axis (a gene the source lacks is a zero column, as a D/E/F file carries every
official gene), 18,400 cells per context with 46 `ntc_id`s × 400, the CURRENT pert panel and manifest.

  bundle X: D = HCT116 (Orion NTC, FiCS; a registry line, strong bank) · E = H1 (Arc 2025 Flex controls, h1-test; absent
            from CCLE) · F = HEK293T (Orion NTC; absent from CCLE, 590's false-HAP1 case)
  bundle Y: D = A549 (Myllia NTC + gene-desert, 10x 5' HT; registry, weak banks) · E = THP-1 (Papalexi ECCITE controls,
            10x 5'; registry, weak bank) · F = Calu-3 (Sunshine NTC, 10x 3', SARS-CoV-2 screen; registry, strong bank)
  Bundle Y's sources have 2,386–4,156 control cells: they are resampled with replacement to 18,400 (seeded); the control
  stages read pseudobulks and library-size quantiles, which resampling preserves.

Writes <out>/{X,Y}/{context_D,E,F.h5ad, gene_names.csv, pert_counts.csv, manifest.json} and a build record with shas.
Usage: uv run python experiments/595-foreign-control-rehearsal/build595.py <out dir>
"""

from __future__ import annotations

import glob
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp

pa.set_memory_pool(pa.system_memory_pool())  # 595: the default pool kept ~12 GB of parquet list columns resident

ROOT = Path(__file__).resolve().parents[2]
OFFICIAL = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
POS = {g: i for i, g in enumerate(OFFICIAL)}
N_CELLS, N_NTC = 18400, 46
SEED = 595


def log(m: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {m}", flush=True)


def colmap(genes: list[str]) -> np.ndarray:
    """dataset gene index -> official column (-1 = not on the official axis); first occurrence of a symbol wins."""
    seen, out = set(), np.full(len(genes), -1, np.int64)
    for j, g in enumerate(genes):
        i = POS.get(g)
        if i is not None and i not in seen:
            seen.add(i)
            out[j] = i
    return out


def remap_csr(m: sp.csr_matrix, cmap: np.ndarray) -> sp.csr_matrix:
    m = m.tocoo()
    keep = cmap[m.col] >= 0
    out = sp.csr_matrix((m.data[keep].astype(np.float32), (m.row[keep], cmap[m.col[keep]])),
                        shape=(m.shape[0], len(OFFICIAL)))
    out.sum_duplicates()
    out.eliminate_zeros()
    out.sort_indices()
    return out


def rows_from_h5ad(path: Path, rows: np.ndarray, var_key: str | None = None) -> tuple[sp.csr_matrix, list[str]]:
    """Selected rows of X (CSR or CSC, h5py; never the whole matrix) and the var names."""
    rows = np.sort(np.asarray(rows))
    with h5py.File(path, "r") as h:
        v = h["var"]
        key = var_key or v.attrs.get("_index", "_index")
        o = v[key]
        if isinstance(o, h5py.Group) and "categories" in o:  # categorical index
            cats = [x.decode() if isinstance(x, bytes) else str(x) for x in o["categories"][:]]
            names = [cats[c] for c in o["codes"][:]]
        elif isinstance(o, h5py.Group) and "values" in o:  # anndata nullable-string-array
            names = [x.decode() if isinstance(x, bytes) else str(x) for x in o["values"][:]]
        else:
            names = [x.decode() if isinstance(x, bytes) else str(x) for x in o[:]]
        X = h["X"]
        enc = X.attrs["encoding-type"]
        n_obs, n_var = (int(s) for s in X.attrs["shape"])
        indptr = X["indptr"][:].astype(np.int64)
        if enc == "csr_matrix":
            parts = []
            for r in rows:
                a, b = indptr[r], indptr[r + 1]
                parts.append((X["data"][a:b], X["indices"][a:b]))
            data = np.concatenate([p[0] for p in parts]).astype(np.float32)
            idx = np.concatenate([p[1] for p in parts]).astype(np.int64)
            ptr = np.concatenate([[0], np.cumsum([len(p[0]) for p in parts])])
            return sp.csr_matrix((data, idx, ptr), shape=(len(rows), n_var)), names
        # CSC: walk column chunks, keep entries whose row is selected
        sel = np.full(n_obs, -1, np.int64)
        sel[rows] = np.arange(len(rows))
        r_all, c_all, d_all = [], [], []
        c0 = 0
        while c0 < n_var:
            c1 = min(n_var, int(np.searchsorted(indptr, indptr[c0] + 20_000_000, side="right")) - 1)
            c1 = max(c1, c0 + 1)
            a, b = indptr[c0], indptr[c1]
            rr = X["indices"][a:b].astype(np.int64)
            dd = X["data"][a:b]
            cc = np.repeat(np.arange(c0, c1), np.diff(indptr[c0:c1 + 1]))
            m = sel[rr] >= 0
            r_all.append(sel[rr[m]]); c_all.append(cc[m]); d_all.append(dd[m].astype(np.float32))
            c0 = c1
        m = sp.coo_matrix((np.concatenate(d_all), (np.concatenate(r_all), np.concatenate(c_all))), shape=(len(rows), n_var))
        return m.tocsr(), names


def obs_col(path: Path, key: str) -> np.ndarray:
    with h5py.File(path, "r") as h:
        o = h["obs"][key]
        if isinstance(o, h5py.Group):
            cats = np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in o["categories"][:]] + ["__NA__"], dtype=object)
            codes = o["codes"][:].astype(np.int64)
            codes[codes < 0] = len(cats) - 1
            return cats[codes]
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in o[:]], dtype=object)


def orion(line: str) -> tuple[sp.csr_matrix, dict]:
    meta = pq.read_table(ROOT / "data/external/orion/metadata/gene_metadata.parquet").to_pandas()
    tok2name = dict(zip(meta.gene_token_id.astype(int), meta.gene_name.astype(str)))
    ntok = int(meta.gene_token_id.max()) + 1
    cmap = np.full(ntok, -1, np.int64)
    seen = set()
    for t in range(ntok):
        i = POS.get(tok2name.get(t, ""))
        if i is not None and i not in seen:
            seen.add(i)
            cmap[t] = i
    files = sorted(glob.glob(str(ROOT / f"data/external/orion/data/{line}_Batch*.parquet")),
                   key=lambda f: int(Path(f).stem.split("Batch")[1]))
    blocks, n, used, integer = [], 0, [], True
    for f in files:
        got = 0
        for rb in pq.ParquetFile(f).iter_batches(batch_size=1024, columns=["gene_token_id", "gene_expression", "gene_target"]):
            keep = np.asarray(rb.column("gene_target").to_pylist(), dtype=object) == "Non-Targeting"
            if not keep.any():
                continue
            rb = rb.filter(pa.array(keep))
            tok = rb.column("gene_token_id")
            val = rb.column("gene_expression")
            offs = tok.offsets.to_numpy()
            toks = tok.values.to_numpy().astype(np.int64)[offs[0]:offs[-1]]
            vals = val.values.to_numpy()[val.offsets.to_numpy()[0]:val.offsets.to_numpy()[-1]]
            offs = offs - offs[0]
            integer &= bool(np.all(vals == np.round(vals)))
            cols = cmap[toks]
            rows = np.repeat(np.arange(rb.num_rows), np.diff(offs))
            k = cols >= 0
            m = sp.csr_matrix((vals[k].astype(np.float32), (rows[k], cols[k])), shape=(rb.num_rows, len(OFFICIAL)))
            m.sum_duplicates()
            blocks.append(m)
            got += rb.num_rows
        n += got
        used.append(Path(f).name)
        log(f"  orion {line}: {Path(f).name} +{got} NTC (total {n})")
        if n >= N_CELLS:
            break
    allm = sp.vstack(blocks).tocsr()
    rng = np.random.default_rng(SEED)
    pick = np.sort(rng.choice(allm.shape[0], N_CELLS, replace=False))
    return allm[pick], {"source": f"X-Atlas/Orion {line}, gene_target == 'Non-Targeting'", "files": used,
                        "pool": int(allm.shape[0]), "resampled": False, "integer_counts": integer}


def h1() -> tuple[sp.csr_matrix, dict]:
    p = ROOT / "data/harness/h1-test.h5ad"
    tg = obs_col(p, "target_gene")
    rows = np.flatnonzero(tg == "non-targeting")
    m, names = rows_from_h5ad(p, rows)
    return remap_csr(m, colmap(names)), {"source": "Arc VCC 2025 H1 (10x Flex), data/harness/h1-test.h5ad controls",
                                         "pool": int(len(rows)), "resampled": False}


def resampled(m: sp.csr_matrix, src: dict) -> tuple[sp.csr_matrix, dict]:
    rng = np.random.default_rng(SEED)
    pick = rng.choice(m.shape[0], N_CELLS, replace=True)
    return m[pick], {**src, "pool": int(m.shape[0]), "resampled": True}


def myllia() -> tuple[sp.csr_matrix, dict]:
    parts = []
    for ch in (1, 2):
        p = ROOT / f"data/external/myllia-a549/GSE337804_A549_ch_{ch}.RDS.h5ad"
        sg = obs_col(p, "sgrna_symbol")
        rows = np.flatnonzero(np.isin(sg, ["non-targeting", "gene_desert"]))
        m, names = rows_from_h5ad(p, rows)
        parts.append(remap_csr(m, colmap(names)))
    return resampled(sp.vstack(parts).tocsr(), {"source": "Myllia A549 GSE337804 (10x 5' HT), non-targeting + gene_desert"})


def scp(fname: str, label: str, src: str) -> tuple[sp.csr_matrix, dict]:
    p = ROOT / "data/external/scperturb" / fname
    pert = obs_col(p, "perturbation")
    rows = np.flatnonzero(pert == "control")
    m, names = rows_from_h5ad(p, rows)
    return resampled(remap_csr(m, colmap(names)), {"source": f"{src} ({label}), perturbation == 'control'"})


def write_context(m: sp.csr_matrix, ctx: str, path: Path) -> dict:
    m = m.tocsr().astype(np.float32)
    m.eliminate_zeros()
    m.sort_indices()
    ntc = [f"non-targeting-{1 + (i // 400)}" for i in range(N_CELLS)]
    obs = pd.DataFrame({"target_gene": pd.Categorical(["non-targeting"] * N_CELLS),
                        "context": pd.Categorical([ctx] * N_CELLS),
                        "ntc_id": pd.Categorical(ntc)}, index=[f"{ctx}_ctrl_{i:06d}" for i in range(N_CELLS)])
    a = ad.AnnData(X=m, obs=obs, var=pd.DataFrame(index=pd.Index(OFFICIAL)))
    a.write_h5ad(path)
    lib = np.asarray(m.sum(1)).ravel()
    q25, q50, q75 = np.percentile(lib, [25, 50, 75])
    return {"cells": int(m.shape[0]), "nnz": int(m.nnz), "genes_with_counts": int((np.asarray(m.sum(0)).ravel() > 0).sum()),
            "lib_q25": float(q25), "lib_median": float(q50), "lib_q75": float(q75), "q75_q25": float(q75 / max(q25, 1)),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def main() -> None:
    out = Path(sys.argv[1])
    plan = {"X": {"D": ("HCT116", lambda: orion("HCT116")), "E": ("H1", h1), "F": ("HEK293T", lambda: orion("HEK293T"))},
            "Y": {"D": ("A549", myllia),
                  "E": ("THP-1", lambda: scp("PapalexiSatija2021_eccite_RNA.h5ad", "THP-1", "Papalexi 2021 ECCITE-seq")),
                  "F": ("Calu-3", lambda: scp("SunshineHein2023.h5ad", "Calu-3", "Sunshine 2023"))}}
    record = {"seed": SEED, "official_axis": "data/challenge/gene_names.csv", "panel": "data/challenge/pert_counts.csv",
              "bundles": {}}
    template = json.loads((ROOT / "experiments/180-def-readiness/work/def-bundle/manifest.json").read_text())
    for b, ctxs in plan.items():
        d = out / b
        d.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / "data/challenge/gene_names.csv", d / "gene_names.csv")
        shutil.copy(ROOT / "data/challenge/pert_counts.csv", d / "pert_counts.csv")
        man = dict(template)
        man["partition"] = "synthetic-def"
        man["panel_id"] = f"vcc2026-mock-foreign-595{b}"
        (d / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
        record["bundles"][b] = {}
        for ctx, (line, fn) in ctxs.items():
            t0 = time.time()
            log(f"bundle {b} context {ctx} = {line}")
            m, src = fn()
            rec = write_context(m, ctx, d / f"context_{ctx}.h5ad")
            rec.update(line=line, **src, seconds=round(time.time() - t0, 1))
            record["bundles"][b][ctx] = rec
            log(f"  wrote {d / f'context_{ctx}.h5ad'}: {rec['cells']} cells, nnz {rec['nnz']}, q75/q25 {rec['q75_q25']:.3f}, "
                f"median lib {rec['lib_median']:.0f}")
            del m
    (out / "build595.json").write_text(json.dumps(record, indent=1) + "\n")
    log(f"wrote {out / 'build595.json'}")


if __name__ == "__main__":
    main()
