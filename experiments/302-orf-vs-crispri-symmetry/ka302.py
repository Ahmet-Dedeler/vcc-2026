"""302 — KA-C: the H1 bank's rows reproduce 279's `Ccross` through the same loader.

Uses NO reserve data (the atlas is not opened).  Recomputes every cross-line pair cosine that
involves `gse295214_h1` from the on-disk banks with 279's own `load_unit`, and compares to
`held279.npz`'s stored `Ccross` rows.  Required to read 0.0e+00.

    PYTHONDONTWRITEBYTECODE=1 uv run python ka302.py
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
WORK = ROOT / "experiments/302-orf-vs-crispri-symmetry/work"
sys.path.insert(0, str(ROOT / "experiments/276-conservation-map"))
sys.path.insert(0, str(ROOT / "experiments/279-znf-state-dependence"))
import lib276 as L  # noqa: E402
import build279 as B  # noqa: E402
import inv279 as I  # noqa: E402
SRC = I.CANDIDATES


def main() -> None:
    z = np.load(ROOT / "experiments/279-znf-state-dependence/work/held279.npz", allow_pickle=True)
    keys = np.asarray(z["keys"]).astype(str).tolist()
    lines = np.asarray(z["lines"]).astype(str).tolist()
    axis = np.asarray(z["axis"]).astype(str).tolist()
    genes = np.asarray(z["genes"]).astype(str).tolist()
    pairs = np.asarray(z["pairs_cross"])
    Ccross = np.asarray(z["Ccross"])
    h1 = keys.index("gse295214_h1")
    sel = [p for p in range(len(pairs)) if h1 in tuple(pairs[p])]
    print(f"axis {len(axis)}, held genes {len(genes)}, cross pairs {len(pairs)}; "
          f"pairs involving gse295214_h1: {len(sel)}", flush=True)

    idxH, UH, _ = B.load_unit(SRC["gse295214_h1"][0], axis, genes)
    print(f"gse295214_h1 unit rows {UH.shape}", flush=True)

    worst = 0.0
    n_cmp = 0
    per_pair = {}
    for p in sel:
        i, j = int(pairs[p][0]), int(pairs[p][1])
        other = keys[j] if i == h1 else keys[i]
        t0 = time.time()
        idxO, UO, _ = B.load_unit(SRC[other][0], axis, genes)
        if i == h1:
            g, c = B.pair_cos(idxH, UH, idxO, UO)
        else:
            g, c = B.pair_cos(idxO, UO, idxH, UH)
        ref = Ccross[p]
        got = np.full(ref.shape, np.nan, np.float32)
        got[g] = c
        both = np.isfinite(ref) | np.isfinite(got)
        d = float(np.nanmax(np.abs(np.nan_to_num(ref[both], nan=0.0)
                                   - np.nan_to_num(got[both], nan=0.0)))) if both.any() else 0.0
        same_n = int(np.isfinite(ref).sum()) == int(np.isfinite(got).sum()) == g.size
        worst = max(worst, d)
        n_cmp += int(g.size)
        per_pair[f"{keys[i]}|{keys[j]}"] = dict(n=int(g.size), max_abs_diff=d, same_n=same_n)
        print(f"  {other:22s} ({lines[j] if i == h1 else lines[i]:14s}) n {g.size:6d}  "
              f"max|diff| {d:.1e}  same_n {same_n}  ({time.time() - t0:.1f}s)", flush=True)
        del idxO, UO

    ok = bool(worst == 0.0 and all(v["same_n"] for v in per_pair.values()))
    print(f"\nKA-C  {len(sel)} pairs, {n_cmp} (pair, gene) cosines, max|diff| vs held279.npz "
          f"Ccross {worst:.1e} -> {'PASS' if ok else 'FAIL'}")
    out = dict(n_pairs=len(sel), n_cosines=n_cmp, max_abs_diff=worst, pass_=ok, per_pair=per_pair,
               rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30)
    (WORK / "ka302_C.json").write_text(json.dumps(out, indent=1))
    print(f"peak RSS {out['rss_gib']:.2f} GiB")


if __name__ == "__main__":
    sys.exit(main())
