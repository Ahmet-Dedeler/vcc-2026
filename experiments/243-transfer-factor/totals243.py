"""243 — stream per-(context, target, gene) panel totals off a shipped h5ad.

214's `block_totals` (experiments/214-cheap-positives/chal214.py lines 48-77), imported by
value: identical loop, identical output keys, plus 196-style per-block cell statistics.
Run ONCE per artifact, strictly serially, under `/usr/bin/time -l`.

KNOWN-ANSWER ROW: `--ka experiments/214-cheap-positives/work/totals-cd4f034.npz` compares the
freshly streamed totals against 214's cached ones and must match on all 900 blocks exactly.

    uv run python totals243.py <arm> <h5ad> [--ka <npz>]
"""
from __future__ import annotations

import hashlib
import sys
import time

import h5py
import numpy as np

from lib243 import ROOT, WORK, rss_gb

N = 400


def block_totals(path):
    """Per-gene panel totals for all 900 blocks, plus nnz and max cell total (214 verbatim,
    with per-block cell-total min/max/mean/sd and per-block nnz added)."""
    with h5py.File(path, "r") as h:
        o, c = h["obs/target_gene"], h["obs/context"]
        tg = np.asarray([x.decode() for x in o["categories"][:]])[o["codes"][:]]
        cx = np.asarray([x.decode() for x in c["categories"][:]])[c["codes"][:]]
        node = h["var/_index"]
        raw = node["values"][:] if isinstance(node, h5py.Group) else node[:]
        genes = np.asarray([g.decode() if isinstance(g, bytes) else str(g) for g in raw])
        G = genes.size
        labels = np.char.add(np.char.add(cx, "|"), tg)
        ptr = h["X/indptr"][:]
        nrow = labels.size
        assert nrow % N == 0
        blocks, tot, nnz, maxcell, cellstat, nnzb = [], [], 0, 0, [], []
        dsd, dsi = h["X/data"], h["X/indices"]
        t0 = time.time()
        for s in range(0, nrow, N):
            lab = labels[s]
            assert (labels[s:s + N] == lab).all(), f"block at {s} is not homogeneous"
            a, b = int(ptr[s]), int(ptr[s + N])
            d = dsd[a:b]
            idx = dsi[a:b]
            t = np.bincount(idx, weights=d, minlength=G).astype(np.int64)
            per = np.diff(ptr[s:s + N + 1])
            rows = np.bincount(np.repeat(np.arange(N), per), weights=d, minlength=N)
            blocks.append(lab)
            tot.append(t)
            nnz += int(b - a)
            nnzb.append(int(b - a))
            maxcell = max(maxcell, int(rows.max()))
            cellstat.append((rows.min(), rows.max(), rows.mean(), rows.std()))
            if (s // N + 1) % 150 == 0:
                print(f"  block {s//N+1}/{nrow//N} [{time.time()-t0:.0f}s] rss {rss_gb():.2f} GB",
                      flush=True)
    return dict(blocks=np.asarray(blocks), totals=np.asarray(tot), genes=genes,
                nnz=nnz, max_cell=maxcell, cellstat=np.asarray(cellstat),
                nnz_block=np.asarray(nnzb))


def main() -> None:
    arm, path = sys.argv[1], ROOT / sys.argv[2]
    ka = None
    if "--ka" in sys.argv:
        ka = ROOT / sys.argv[sys.argv.index("--ka") + 1]
    t0 = time.time()
    r = block_totals(path)
    out = WORK / f"totals-{arm}.npz"
    np.savez_compressed(out, **r)
    sha = hashlib.sha256(r["totals"].tobytes()).hexdigest()
    print(f"{arm}: {path.name}  blocks {r['blocks'].size}  genes {r['genes'].size}  "
          f"nnz {r['nnz']}  max_cell {r['max_cell']}  totals sha256 {sha[:16]}  "
          f"[{time.time()-t0:.0f}s]  peak rss {rss_gb():.2f} GB")
    if ka is not None:
        z = np.load(ka, allow_pickle=True)
        same_blocks = bool(np.array_equal(np.asarray(z["blocks"]).astype(str), r["blocks"]))
        same = int(sum(np.array_equal(a, b) for a, b in zip(z["totals"], r["totals"])))
        nnz_ok = int(z["nnz"]) == r["nnz"]
        print(f"KA vs {ka.relative_to(ROOT)}: blocks identical {same_blocks}; totals identical on "
              f"{same}/{r['blocks'].size} blocks; nnz {int(z['nnz'])} vs {r['nnz']} "
              f"{'PASS' if (same_blocks and same == r['blocks'].size and nnz_ok) else 'FAIL'}")
        if not (same_blocks and same == r["blocks"].size and nnz_ok):
            raise SystemExit(1)


if __name__ == "__main__":
    main()
