"""Step 9: seed repeats on the leading arms.

The sweep separates the best per-gene policy from the plain global one by 0.0006 on avg6,
which is not obviously bigger than the emitter's own seed noise. Three extra independent
emission seeds per arm say whether that gap is real. It decides the recommendation: if the
arms are indistinguishable, ship the one with no `keep_top` at all.
"""
from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib093 as L  # noqa: E402
import s4_sweep as S  # noqa: E402

OUT = os.path.join(L.SCRATCH, "seeds093")
SEEDS = (11, 12, 13)
ARMS = [("global_s1.0", 0, 1.0, "eff", L.PROD_SPREAD),
        ("global_s0.85", 0, 0.85, "eff", L.PROD_SPREAD),
        ("effK600_s1.0", 600, 1.0, "eff", L.PROD_SPREAD),
        ("effK850_s1.0", 850, 1.0, "eff", L.PROD_SPREAD),
        ("tightK850_s1.0", 850, 1.0, "eff", 0.0),
        ("padjK850_s1.0", 850, 1.0, "padj", L.PROD_SPREAD)]
CHUNK, WORKERS = 4, 5


def _work(chunk):
    gate = S._G
    out = {}
    for name, *_ in ARMS:
        for sd in SEEDS:
            out[f"{name}__s{sd}"] = tuple(np.empty((len(chunk), gate.G)) for _ in range(3))
    for ci, i in enumerate(chunk):
        for name, k, s, key, keep in ARMS:
            vec = S._spread_vector(i, k, s, key, keep)
            for sd in SEEDS:
                rng = np.random.default_rng(sd * 1_000_003 + i)
                block = L.emit_block(S._TOT[i], vec, rng)
                if not np.array_equal(block.sum(axis=0).astype(np.int64), S._TOT[i]):
                    raise AssertionError(f"panel totals moved: {i} {name} seed {sd}")
                p, q, lf = gate.de(block)
                del block
                tag = f"{name}__s{sd}"
                out[tag][0][ci] = p
                out[tag][1][ci] = q
                out[tag][2][ci] = lf
    return chunk, out


def main():
    os.makedirs(OUT, exist_ok=True)
    gate = L.Gate()
    gate.ctrl()
    n_p = gate.perts.size
    tags = [f"{n}__s{sd}" for n, *_ in ARMS for sd in SEEDS]
    acc = {t: [np.empty((n_p, gate.G)) for _ in range(3)] for t in tags}
    chunks = [np.arange(a, min(a + CHUNK, n_p)) for a in range(0, n_p, CHUNK)]
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(WORKERS, initializer=S._init) as ex:
        for chunk, res in ex.map(_work, chunks):
            for t in res:
                for f in range(3):
                    acc[t][f][chunk] = res[t][f]
            done += len(chunk)
            print(f"  {done}/{n_p}  {time.time() - t0:.0f}s", flush=True)
    for t in acc:
        np.savez(os.path.join(OUT, f"{t}.npz"), pred_p=acc[t][0],
                 pred_padj=acc[t][1], pred_lfc=acc[t][2])
    print(f"done in {time.time() - t0:.0f}s -> {OUT}")


if __name__ == "__main__":
    main()
