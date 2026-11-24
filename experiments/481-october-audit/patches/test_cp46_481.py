"""481 F2 known-answer and sensitivity rows for the proposed Cp46 size normalisation (`cp46_fix481.py`).

Rows (context D = A's controls, 46 NTC constructs x 400 cells):
  KA      as shipped: the fixed Cp46 must equal ntclevel417's CP46_ABC["A"] = 14.439751080335085 bit for bit;
  pairs   constructs merged in pairs (23 x 800 cells): current rule vs fixed;
  halves  each construct split in two (92 x 200 cells): current vs fixed;
  uneven  the 18,400 cells re-assigned to 46 constructs of uneven size (seeded, 200..700 cells): current vs fixed.
The "current rule" is the fixed function with cells_per_pert set to each construct's own size (factor 1), i.e.
ntclevel417's arithmetic. Peak ~1.8 GB.
usage: uv run python patches/test_cp46_481.py -> work/test_cp46_481.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csr_matrix
from cell_eval2.moments import jackknife_correction

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cp46_fix481 import cp46_rows  # noqa: E402

EXP = HERE.parent
CP46_ABC_A = 14.439751080335085
path = EXP.parent / "180-def-readiness/work/def-bundle/context_D.h5ad"
with h5py.File(path, "r") as h:
    X = csr_matrix((h["X"]["data"][:], h["X"]["indices"][:], h["X"]["indptr"][:]), shape=tuple(h["X"].attrs["shape"]))
    g = h["obs"]["ntc_id"]
    ids = np.asarray(g["categories"][:]).astype(str)[np.asarray(g["codes"][:])]


def jk(block):
    return float(jackknife_correction(csr_matrix(block), np.zeros(block.shape[0], np.intp), 1, 50000.0, chunk=64)[0])


def current(ids_):
    vals = [jk(X[np.flatnonzero(ids_ == n)]) for n in dict.fromkeys(ids_)]
    return float(np.mean(vals))


out = {}
m, sd, k, sizes = cp46_rows(X, ids, 400, jk)
out["KA_as_shipped"] = {"fixed": m, "ntclevel417_CP46_ABC_A": CP46_ABC_A, "bit_identical": m == CP46_ABC_A,
                        "n_constructs": k, "sizes": sorted(set(sizes))}
names = list(dict.fromkeys(ids))
pairs = ids.copy()
for j in range(0, len(names), 2):
    pairs[ids == names[j + 1]] = names[j]
halves = ids.copy().astype(object)
for n in names:
    r = np.flatnonzero(ids == n)
    halves[r[200:]] = n + "_b"
halves = halves.astype(str)
rng = np.random.default_rng(481)
w = rng.uniform(200, 700, 46)
w = np.floor(w / w.sum() * ids.size).astype(int)
w[: ids.size - w.sum()] += 1
uneven = np.repeat(np.array([f"u{i}" for i in range(46)]), w)[rng.permutation(ids.size)]
for tag, lab in (("pairs_800", pairs), ("halves_200", halves), ("uneven_200_700", uneven)):
    cur = current(lab)
    fx = cp46_rows(X, lab, 400, jk)[0]
    out[tag] = {"current_rule": cur, "fixed": fx, "current_over_true": cur / CP46_ABC_A, "fixed_over_true": fx / CP46_ABC_A,
                "cliff_scale_error_current": float(np.sqrt(cur / CP46_ABC_A)) - 1.0,
                "cliff_scale_error_fixed": float(np.sqrt(fx / CP46_ABC_A)) - 1.0}
    print(tag, {k2: round(v, 4) for k2, v in out[tag].items()}, flush=True)
print("KA", out["KA_as_shipped"], flush=True)
(EXP / "work/test_cp46_481.json").write_text(json.dumps(out, indent=1))
if not out["KA_as_shipped"]["bit_identical"]:
    raise SystemExit("KA FAIL")
print("KA PASS", flush=True)
