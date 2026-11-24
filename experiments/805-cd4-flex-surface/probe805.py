"""805 probe (runs on EC2 only): the HDF5 layout of the GWCD4i cell-level files, read by range requests from the public
bucket (anonymous, us-east-1). Reads ONLY structure, obs columns, var and X/indptr, plus a timed read of 64 random rows'
X/data + X/indices (the values are discarded, never summarised). No DE, no expression statistic of any kind.

    uv run --with s3fs python experiments/805-cd4-flex-surface/probe805.py <out.json> D4_Rest [D1_Rest ...]
"""
from __future__ import annotations

import json
import sys
import time

import h5py
import numpy as np
import s3fs

BUCKET = "genome-scale-tcell-perturb-seq/marson2025_data"


def describe(node, depth=0, out=None, prefix=""):
    out = {} if out is None else out
    for k, v in node.items():
        p = f"{prefix}/{k}"
        if isinstance(v, h5py.Dataset):
            out[p] = dict(shape=list(v.shape), dtype=str(v.dtype), chunks=list(v.chunks) if v.chunks else None,
                          compression=v.compression, compression_opts=str(v.compression_opts),
                          attrs={a: str(v.attrs[a])[:200] for a in v.attrs})
        else:
            out[p] = dict(group=True, attrs={a: str(v.attrs[a])[:200] for a in v.attrs})
            if depth < 3:
                describe(v, depth + 1, out, p)
    return out


def main() -> None:
    outp, names = sys.argv[1], sys.argv[2:]
    fs = s3fs.S3FileSystem(anon=True, default_block_size=16 << 20)
    res = {}
    for name in names:
        t0 = time.time()
        key = f"{BUCKET}/{name}.assigned_guide.h5ad"
        f = fs.open(key, "rb", block_size=16 << 20, cache_type="readahead")
        h = h5py.File(f, "r")
        r = dict(size=fs.size(key), layout=describe(h))
        r["t_open"] = time.time() - t0
        t1 = time.time()
        ip = h["X/indptr"][:]
        r["t_indptr"] = time.time() - t1
        r["n_cells"] = int(len(ip) - 1)
        r["nnz"] = int(ip[-1])
        r["nnz_per_cell_median"] = float(np.median(np.diff(ip)))
        t2 = time.time()
        obs = {}
        for col in ("guide_type", "low_quality", "perturbed_gene_name", "lane_id", "guide_id"):
            node = h["obs"].get(col)
            if node is None:
                continue
            if isinstance(node, h5py.Group):
                cats = node["categories"][:]
                codes = node["codes"][:]
                vals, cnt = np.unique(codes, return_counts=True)
                obs[col] = dict(kind="categorical", n_categories=int(len(cats)),
                                top=[(cats[v].decode() if isinstance(cats[v], bytes) else str(cats[v]), int(c))
                                     for v, c in sorted(zip(vals, cnt), key=lambda z: -z[1])[:8] if v >= 0])
            else:
                x = node[:]
                if x.dtype == bool or x.dtype.kind in "iu":
                    vals, cnt = np.unique(x, return_counts=True)
                    obs[col] = dict(kind=str(x.dtype), counts={str(a): int(b) for a, b in zip(vals[:10], cnt[:10])})
                else:
                    obs[col] = dict(kind=str(x.dtype), first=[str(y) for y in x[:3]])
        r["obs_summary"] = obs
        r["t_obs"] = time.time() - t2
        gn = h["var"].get("gene_name")
        if gn is not None:
            if isinstance(gn, h5py.Group):
                g = gn["categories"][:][gn["codes"][:]]
            else:
                g = gn[:]
            r["var_gene_name_first"] = [x.decode() if isinstance(x, bytes) else str(x) for x in g[:5]]
            r["n_var"] = int(len(g))
        rng = np.random.default_rng(805)
        rows = np.sort(rng.choice(len(ip) - 1, size=64, replace=False))
        t3 = time.time()
        for i in rows:
            a, b = int(ip[i]), int(ip[i + 1])
            _ = h["X/data"][a:b]
            _ = h["X/indices"][a:b]
        r["t_64_random_rows"] = time.time() - t3
        res[name] = r
        h.close()
        f.close()
        print(json.dumps({name: {k: v for k, v in r.items() if k != "layout"}}, default=str), flush=True)
    json.dump(res, open(outp, "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
