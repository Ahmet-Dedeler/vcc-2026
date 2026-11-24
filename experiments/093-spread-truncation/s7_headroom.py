"""Step 7: how much headroom is there PAST the `spread` clip?

`emit_counts` clips `spread` to 1.0, i.e. the scattered part is drawn Poisson. That is the
noisiest a total-preserving scatter can be *within the current emitter*, and the sweep shows
it only takes median `n_pred` from 5360 to ~2440 -- well short of the ~850 the 092
counterfactual wants. The obvious continuation is an OVER-Poisson scatter: draw the random
part gamma-Poisson at dispersion `phi` and repair the gene total the same exact way, so the
pseudobulk invariance -- the whole reason this channel is safe -- is untouched.

This probe is deliberately OUTSIDE `src/`: it measures whether extending the emitter is
worth a follow-up, and does not change production semantics. `phi = 0` reproduces the
`spread = 1.0` arm exactly.
"""
from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib093 as L  # noqa: E402

OUT = os.path.join(L.SCRATCH, "sweep093nb")
PHIS = (0.0, 0.05, 0.15, 0.4, 1.0)
CHUNK = 4
WORKERS = 5


def emit_nb(totals: np.ndarray, phi: float, rng, n: int = L.N_EMIT) -> np.ndarray:
    """The quasi emitter at spread 1.0 with a gamma-Poisson scatter instead of Poisson.

    Structure copied from `donor_transfer.emit_counts`'s quasi branch, including the exact
    residual repair, so the ONLY difference from `spread = 1.0` is the scatter's variance:
    `mu + phi * mu^2` against Poisson's `mu`.
    """
    totals = np.asarray(totals, dtype=np.int64)
    G = totals.size
    block = np.zeros((n, G), dtype=np.int64)
    rz = np.flatnonzero(totals > 0)
    mu = np.broadcast_to(totals[rz] / n, (n, rz.size))
    if phi > 0:
        draw = rng.poisson(rng.gamma(shape=1.0 / phi, scale=mu * phi))
    else:
        draw = rng.poisson(mu)
    resid = totals[rz] - draw.sum(0)
    for j in np.flatnonzero(resid):
        d = int(resid[j])
        if d > 0:
            np.add.at(draw[:, j], rng.integers(0, n, d), 1)
        else:
            for _ in range(-d):
                nz = np.flatnonzero(draw[:, j] > 0)
                if nz.size == 0:
                    break
                draw[rng.choice(nz), j] -= 1
    block[:, rz] = draw
    return block.astype(np.float64)


_G = None
_TOT = None


def _init():
    global _G, _TOT
    _G = L.Gate()
    _G.ctrl()
    _TOT = np.load(os.path.join(L.SCRATCH, "totals093.npy"))


def _work(chunk):
    out = {f"nb_phi{p}": tuple(np.empty((len(chunk), _G.G)) for _ in range(3))
           for p in PHIS}
    for ci, i in enumerate(chunk):
        for p in PHIS:
            rng = np.random.default_rng(70_000 + i)
            block = emit_nb(_TOT[i], p, rng)
            if not np.array_equal(block.sum(axis=0).astype(np.int64), _TOT[i]):
                raise AssertionError(f"panel totals moved: pert {i} phi {p}")
            pv, q, lf = _G.de(block)
            del block
            out[f"nb_phi{p}"][0][ci] = pv
            out[f"nb_phi{p}"][1][ci] = q
            out[f"nb_phi{p}"][2][ci] = lf
    return chunk, out


def main():
    os.makedirs(OUT, exist_ok=True)
    gate = L.Gate()
    gate.ctrl()
    n_p = gate.perts.size
    acc = {f"nb_phi{p}": [np.empty((n_p, gate.G)) for _ in range(3)] for p in PHIS}
    chunks = [np.arange(a, min(a + CHUNK, n_p)) for a in range(0, n_p, CHUNK)]
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(WORKERS, initializer=_init) as ex:
        for chunk, res in ex.map(_work, chunks):
            for name in res:
                for f in range(3):
                    acc[name][f][chunk] = res[name][f]
            done += len(chunk)
            print(f"  {done}/{n_p}  {time.time() - t0:.0f}s", flush=True)
    for name in acc:
        np.savez(os.path.join(OUT, f"{name}.npz"), pred_p=acc[name][0],
                 pred_padj=acc[name][1], pred_lfc=acc[name][2])
    print(f"done in {time.time() - t0:.0f}s -> {OUT}")


if __name__ == "__main__":
    main()
