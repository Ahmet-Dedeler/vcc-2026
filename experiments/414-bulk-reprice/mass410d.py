"""414: the SYVN1 context-C block in 410c vs 410d: mass share of the top genes (is the energy blow-up a render mass explosion?)"""
import os, sys, json
from pathlib import Path
HERE = Path(__file__).resolve().parent
os.environ.setdefault("NUMBA_CACHE_DIR", str(HERE / "work" / "numba-cache"))
import numpy as np, h5py
sys.path.insert(0, str(HERE))
from decomp414 import str_col
ROOT = HERE.parents[1]
out = {}
for tag, p in (("410c", "build_c/the-kettle-applied-for-tenure.h5ad"), ("410d", "build_d/the-lighthouse-keeper-went-remote.h5ad")):
    with h5py.File(ROOT / "experiments/410-pertscale-revisit/work" / p, "r") as h:
        genes = str_col(h, "var/_index"); ctx = str_col(h, "obs/context"); tg = str_col(h, "obs/target_gene")
        indptr = np.asarray(h["X"]["indptr"][:])
        for c in "ABC":
            rows = np.flatnonzero((ctx == c) & (tg == "SYVN1")); lo, hi = int(indptr[rows[0]]), int(indptr[rows[-1] + 1])
            tot = np.bincount(h["X"]["indices"][lo:hi].astype(np.int64), weights=h["X"]["data"][lo:hi].astype(np.float64), minlength=genes.size)
            top = np.argsort(-tot)[:4]
            out[f"{tag}_{c}"] = [(str(genes[g]), round(float(tot[g] / tot.sum()), 4)) for g in top]
            print(tag, c, "SYVN1 block: top-4 genes' share of block counts", out[f"{tag}_{c}"], flush=True)
(HERE / "work" / "mass410d.json").write_text(json.dumps(out, indent=1))
