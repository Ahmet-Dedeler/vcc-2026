"""481: how ntclevel417's Cp46 (the re-level's claim) depends on the NTC construct size.

ntclevel417 takes the mean over the context's NTC constructs of cell_eval2's delete-1 jackknife on each construct's cells,
as the proxy for the reference's per-perturbation sampling noise C_p (400-cell groups). The preflight checks the NUMBER of
NTC ids against the manifest, never the cells per construct. This recomputes Cp on context D (= A) with the constructs
as shipped (400 cells), halved (first 200 cells of each), and merged in pairs (800 cells).
usage: uv run python cpsize481.py [ctx] -> work/cpsize481.json   (peak ~1.5 GB)
"""
import json, sys, time
from pathlib import Path
import h5py, numpy as np
from scipy.sparse import csr_matrix
from cell_eval2.moments import jackknife_correction

HERE = Path(__file__).resolve().parent
ctx = sys.argv[1] if len(sys.argv) > 1 else "D"
path = HERE.parents[0] / "180-def-readiness/work/def-bundle" / f"context_{ctx}.h5ad"
with h5py.File(path, "r") as h:
    X = csr_matrix((h["X"]["data"][:], h["X"]["indices"][:], h["X"]["indptr"][:]), shape=tuple(h["X"].attrs["shape"]))
    g = h["obs"]["ntc_id"]
    ids = np.asarray(g["categories"][:]).astype(str)[np.asarray(g["codes"][:])]
TS = 50000.0
def cp(rows):
    return float(jackknife_correction(X[rows], np.zeros(len(rows), np.intp), 1, TS, chunk=64)[0])
names = list(dict.fromkeys(ids))
rows = {n: np.flatnonzero(ids == n) for n in names}
t0 = time.time()
c400 = [cp(rows[n]) for n in names]
c200 = [cp(rows[n][:200]) for n in names]
c800 = [cp(np.concatenate([rows[names[2 * k]], rows[names[2 * k + 1]]])) for k in range(len(names) // 2)]
out = {"ctx": ctx, "n_constructs": len(names), "Cp_400": float(np.mean(c400)), "Cp_200": float(np.mean(c200)),
       "Cp_800": float(np.mean(c800)), "ratio_200_over_400": float(np.mean(c200) / np.mean(c400)),
       "ratio_800_over_400": float(np.mean(c800) / np.mean(c400)), "seconds": round(time.time() - t0, 1)}
print(json.dumps(out))
(HERE / "work" / f"cpsize481_{ctx}.json").write_text(json.dumps(out, indent=1))
