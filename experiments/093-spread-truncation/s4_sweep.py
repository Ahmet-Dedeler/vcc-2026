"""Step 4: the measurement. Re-emit every perturbation under every spread policy and run
the exact DE path on the result.

Policies
--------
global s            every gene at spread `s` (K = 0). Production is the s = 0.25 row, so
                    that row doubles as the reproduction control.
effK<K>_s<S>        top K genes by |predicted effect| held at the production 0.25, the rest
                    at S. The form the brief specifies.
padjK<K>_s<S>       ORACLE arm: the top K by the production run's own adjusted p-value.
                    Not buildable (it needs the emission it is choosing for), but it bounds
                    what ANY per-gene selector could win, which is the number that decides
                    whether a better pre-emission key is worth chasing.
tightK<K>_s<S>      top K by |effect| pushed all the way to the variance floor (0.0)
                    instead of 0.25, the rest at S -- widens the gap the channel can open.

Every arm asserts the emitted panel totals are bit-identical to production's, which is what
makes `pds` and `mse` invariant.
"""
from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib093 as L  # noqa: E402

OUT = os.path.join(L.SCRATCH, "sweep093")
CHUNK = 4
WORKERS = 5


def policies():
    pol = [(f"global_s{s}", 0, s, "eff", L.PROD_SPREAD) for s in
           (0.25, 0.4, 0.55, 0.7, 0.85, 1.0)]
    for s in (0.6, 0.8, 1.0):
        for k in (400, 600, 850, 1200, 2000):
            pol.append((f"effK{k}_s{s}", k, s, "eff", L.PROD_SPREAD))
        pol.append((f"effP15_s{s}", -0.15, s, "eff", L.PROD_SPREAD))
    for k in (850, 2000):
        pol.append((f"padjK{k}_s1.0", k, 1.0, "padj", L.PROD_SPREAD))
        pol.append((f"tightK{k}_s1.0", k, 1.0, "eff", 0.0))
    return pol


_G = None
_TOT = None
_KEYS = None


def _init():
    global _G, _TOT, _KEYS
    _G = L.Gate()
    _G.ctrl()
    _TOT = np.load(os.path.join(L.SCRATCH, "totals093.npy"))
    d = np.load(os.path.join(L.E092, "de_tables.npz"), allow_pickle=True)
    pq = d["pred_padj"]
    pp = d["pred_p"]
    pl = d["pred_lfc"]
    genes = d["genes"].astype(str)
    fr = np.empty(genes.size, np.int64)
    fr[np.argsort(genes)] = np.arange(genes.size)
    _KEYS = dict(pq=pq, pp=pp, pl=pl, fr=fr,
                 gpos={g: i for i, g in enumerate(genes)},
                 perts=d["perts"].astype(str))


def _spread_vector(i: int, k, s: float, key: str, keep_spread: float) -> np.ndarray:
    """Per-gene spread over the FULL gene axis."""
    gate = _G
    vec = np.full(gate.G_full, float(s))
    if k == 0:
        return vec
    if key == "eff":
        score = L.rank_key(_TOT[i], gate)          # |predicted effect|, gated axis
        n_call = None
    else:                                           # oracle: our own adjusted p-value
        pq = _KEYS["pq"][i]
        score = -pq
        n_call = int((pq < L.ALPHA).sum())
    if k < 0:                                       # fraction of our own call list
        pq = _KEYS["pq"][i]
        n_call = int((pq < L.ALPHA).sum()) if n_call is None else n_call
        kk = int(round(-k * n_call))
    else:
        kk = int(k)
    kk = min(kk, score.size)
    top = np.argpartition(-score, kk - 1)[:kk]
    vec[gate.keep_idx[top]] = float(keep_spread)
    return vec


def _work(chunk):
    gate = _G
    pol = policies()
    out = {}
    for name, *_ in pol:
        out[name] = (np.empty((len(chunk), gate.G)), np.empty((len(chunk), gate.G)),
                     np.empty((len(chunk), gate.G)))
    for ci, i in enumerate(chunk):
        for name, k, s, key, keep in pol:
            vec = _spread_vector(i, k, s, key, keep)
            rng = np.random.default_rng(90_000 + i)
            block = L.emit_block(_TOT[i], vec, rng)
            if not np.array_equal(block.sum(axis=0).astype(np.int64), _TOT[i]):
                raise AssertionError(f"panel totals moved: pert {i} policy {name}")
            p, q, lf = gate.de(block)
            del block
            out[name][0][ci] = p
            out[name][1][ci] = q
            out[name][2][ci] = lf
    return chunk, out


def main():
    os.makedirs(OUT, exist_ok=True)
    gate = L.Gate()
    gate.ctrl()          # warm the tie-sum cache before forking workers
    n_p = gate.perts.size
    pol = policies()
    print(f"{len(pol)} policies x {n_p} perturbations, {WORKERS} workers", flush=True)
    acc = {name: [np.empty((n_p, gate.G)) for _ in range(3)] for name, *_ in pol}
    chunks = [np.arange(a, min(a + CHUNK, n_p)) for a in range(0, n_p, CHUNK)]
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(WORKERS, initializer=_init) as ex:
        for chunk, res in ex.map(_work, chunks):
            for name in res:
                for f in range(3):
                    acc[name][f][chunk] = res[name][f]
            done += len(chunk)
            el = time.time() - t0
            print(f"  {done}/{n_p} perts   {el:.0f}s elapsed, "
                  f"eta {el / done * (n_p - done):.0f}s", flush=True)
    for name in acc:
        np.savez(os.path.join(OUT, f"{name}.npz"),
                 pred_p=acc[name][0], pred_padj=acc[name][1], pred_lfc=acc[name][2])
    print(f"done in {time.time() - t0:.0f}s -> {OUT}")


if __name__ == "__main__":
    main()
