#!/usr/bin/env python3
"""658 pre-slot check: do two prediction h5ads carry the SAME bulk, i.e. identical per-(context, target, gene) count sums?
That is exactly the condition for identical pds / mse on the board (both read the log-pseudobulk of summed counts). With
--exact, also require the whole matrix (cell order, indptr, indices, data) and the obs labels to be identical (524dr2's N0
branch: same recipe, applied tree => same content).

  uv run python experiments/658-prereg-1007-builds/bulkid658.py A.h5ad B.h5ad [--exact] [--out J]

Use: ALLSW vs 498p2g's package restored with scripts/restore_prediction_h5ad.py (F1 checked before the slot); 524dr2 vs 524dr
the same way. Streams row blocks from h5py (no AnnData load); sums are int64 per (context, target) x 18,533 genes, ~0.13 GB per
file. Exit 0 = identical, 1 = different, 2 = the two files are not comparable (different gene axis or row labels)."""
import json
import sys
from pathlib import Path

import h5py
import numpy as np
from anndata.io import read_elem

BLOCK = 4000                       # rows per read: ~1.2 GB peak on a full 12 GB submission (20000 rows: 5.4 GB)


def obs_col(f, name):
    """an obs column as str, whatever the encoding (categorical, plain or nullable strings): anndata's element reader."""
    return read_elem(f["obs"])[name].astype(str).to_numpy()


def var_names(f):
    return read_elem(f["var"]).index.astype(str).to_numpy()


def sums(path):
    with h5py.File(path, "r") as f:
        X = f["X"]
        indptr = X["indptr"]
        n, g = len(indptr) - 1, len(var_names(f))
        ctx, tg = obs_col(f, "context"), obs_col(f, "target_gene")
        keys = sorted(set(zip(ctx, tg)))
        kid = {k: i for i, k in enumerate(keys)}
        row_key = np.array([kid[(c, t)] for c, t in zip(ctx, tg)])
        S = np.zeros((len(keys), g), dtype=np.int64)
        ip = np.asarray(indptr[...])
        for a in range(0, n, BLOCK):
            b = min(n, a + BLOCK)
            lo, hi = int(ip[a]), int(ip[b])
            ind = np.asarray(X["indices"][lo:hi])
            dat = np.asarray(X["data"][lo:hi])
            rows = np.repeat(np.arange(a, b), np.diff(ip[a:b + 1]))
            np.add.at(S, (row_key[rows], ind), np.rint(dat).astype(np.int64))
        return keys, S, var_names(f), ctx, tg


def exact(pa, pb):
    with h5py.File(pa, "r") as fa, h5py.File(pb, "r") as fb:
        for col in ("context", "target_gene"):
            if not np.array_equal(obs_col(fa, col), obs_col(fb, col)):
                return f"obs {col} differs"
        ia, ib = np.asarray(fa["X"]["indptr"][...]), np.asarray(fb["X"]["indptr"][...])
        if not np.array_equal(ia, ib):
            return f"indptr differs (first row {int(np.argmax(ia != ib))})"
        n = len(ia) - 1
        for a in range(0, n, BLOCK):
            b = min(n, a + BLOCK)
            lo, hi = int(ia[a]), int(ia[b])
            for k in ("indices", "data"):
                if not np.array_equal(np.asarray(fa["X"][k][lo:hi]), np.asarray(fb["X"][k][lo:hi])):
                    return f"X {k} differs in rows {a}..{b}"
    return None


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
    args = [a for a in args if a != out]
    pa, pb = args[:2]
    ka, Sa, ga, *_ = sums(pa)
    kb, Sb, gb, *_ = sums(pb)
    rec = {"a": pa, "b": pb}
    if not np.array_equal(ga, gb) or ka != kb:
        rec["comparable"] = False
        print(f"NOT COMPARABLE: gene axis equal {np.array_equal(ga, gb)}, (context, target) keys equal {ka == kb}")
        code = 2
    else:
        d = Sa - Sb
        rec.update(comparable=True, cells_differing=int((d != 0).sum()), abs_sum=int(np.abs(d).sum()),
                   blocks_differing=int((d != 0).any(1).sum()), n_blocks=len(ka))
        print(f"bulk: {rec['cells_differing']} (block, gene) sums differ in {rec['blocks_differing']} / {len(ka)} blocks; "
              f"sum |d| = {rec['abs_sum']}")
        code = 0 if rec["cells_differing"] == 0 else 1
        if "--exact" in sys.argv:
            why = exact(pa, pb)
            rec["exact"] = why is None
            rec["exact_reason"] = why
            print("exact content: " + ("IDENTICAL" if why is None else f"DIFFERENT ({why})"))
            code = code if why is None else 1
    print("BULK IDENTICAL" if code == 0 else "DIFFERENT" if code == 1 else "")
    if out:
        Path(out).write_text(json.dumps(rec, indent=1))
    sys.exit(code)


if __name__ == "__main__":
    main()
