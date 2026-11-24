#!/usr/bin/env python3
"""Known-answer / invariant tests for the 423 rank-preserving magnitude pass on real emitted H1 blocks.

Uses the blocks dumped by `run423.py 0 --dump 8` (work/explore/block_s0_0*.npz; base_s0 construction, cap15 bulk).
  1. k = 1.0 is a no-op in every mode (output bit-identical to the input block).
  2. every mode: per-gene block totals and per-cell libraries bit-identical; no negative counts; no lfc sign flips;
     untouched DE-universe genes keep their exact rank sum and tie term.
  3. ex / zs: the rank sum R1 of EVERY DE-universe gene is identical (scorer arithmetic, rp423.r1_tie).
  4. pc: every carrier cell keeps its exact position among the controls for every gene (only the ballast moves).
usage: uv run python test_rp423.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run423 as RN  # noqa: E402  (paths)
import rp423  # noqa: E402

E = RN.E


def load(i):
    z = np.load(HERE / f"work/explore/block_s0_{i:02d}.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    return X, int(z["n_car"])


def main():
    r, cc, _ = E.load_reference()
    fails = 0
    for i in (1, 6):
        X, n_car = load(i)
        s, j, p, q, l = E.block_stats(X, cc)
        R0, T0 = rp423.r1_tie(X, cc)
        for mode in rp423.MODES:
            X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant(f"all_k100_{mode}"))
            ok = np.array_equal(X2, X)
            print(f"[{'PASS' if ok else 'FAIL'}] target {i} k=1.0 {mode}: bit-identical no-op")
            fails += not ok
        for v in ("all_k140", "dn_k250", "all_k140_zs", "all_k140_ex", "dn_k250_ex"):
            X2, dg = rp423.postpass(X, n_car, cc, l, rp423.parse_variant(v))
            touched = dg["_touched"]
            s2, j2, p2, q2, l2 = E.block_stats(X2, cc)
            R1, T1 = rp423.r1_tie(X2, cc)
            fin = np.isfinite(l) & np.isfinite(l2)
            checks = {
                "per-gene totals": np.array_equal(X2.sum(0), X.sum(0)),
                "libraries": np.array_equal(X2.sum(1), X.sum(1)),
                "non-negative": X2.min() >= 0,
                "no sign flip": not np.any(np.sign(l[fin]) != np.sign(l2[fin])),
                "untouched R1+tie exact": np.array_equal(R1[~touched], R0[~touched]) and np.array_equal(T1[~touched], T0[~touched]),
                "something moved": bool(touched.any()),
            }
            if v.endswith(("_ex", "_zs")):
                checks["R1 identical, every gene"] = np.array_equal(R1, R0)
                checks["calls identical"] = np.array_equal(q < .05, q2 < .05)
            else:
                Lr = X.sum(1) / 1e6
                keep = np.flatnonzero(cc.keep)
                same = True
                for jj in np.flatnonzero(touched):
                    a0 = rp423.ranks_of(X[:n_car, keep[jj]], Lr[:n_car], cc.cs[jj])
                    a1 = rp423.ranks_of(X2[:n_car, keep[jj]], Lr[:n_car], cc.cs[jj])
                    same &= np.array_equal(a0, a1)
                checks["carrier positions exact (pc)"] = bool(same)
            for name, okc in checks.items():
                print(f"[{'PASS' if okc else 'FAIL'}] target {i} {v}: {name}")
                fails += not okc
    print(f"{fails} failed")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
