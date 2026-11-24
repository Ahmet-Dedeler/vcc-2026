#!/usr/bin/env python3
"""383 — score STATE H1 predictions through the 323 emitter+scorer, my dir only.

Inputs:  work/state_pred_h1.h5ad (STATE predicted cells, obs[target_gene], log1p X
         in support 18,080-gene order) — built by `state tx infer` on work/templates/h1*.
Steps:   per-target log1p delta -> effect matrix (support gene axis) -> map_axis to
         r[genes] -> 323 emission loop (rc_hybrid_block, production seeds/settings,
         copied call pattern from evaluate323.emit_and_score) -> block_stats -> score().
         Writes work/metrics-state_h1.json + per-target CSV + blend cosines.
Blend:   per-target cosine(STATE, cap15), cosine(STATE, truth), cosine(cap15, truth),
         cosine(50/50 fixed-norm blend, truth) on the scorer keep-universe genes.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("POLARS_MAX_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
for _p in (str(BASE), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np
import polars as pl

import evaluate323 as E  # noqa: E402
import run323 as R  # noqa: E402
from vcc2026 import rc_emit  # noqa: E402
from vcc2026.donor_transfer import emit_counts  # noqa: E402


def state_effect_matrix(pred_h5ad: str) -> tuple[np.ndarray, list, list]:
    import anndata as ad

    a = ad.read_h5ad(pred_h5ad)
    tg = a.obs["target_gene"].astype(str).to_numpy()
    X = a.X.toarray() if not isinstance(a.X, np.ndarray) else a.X
    X = np.asarray(X, dtype=np.float64)
    genes = [str(g) for g in a.var_names]
    perts = sorted(set(tg) - {"non-targeting"})
    cm = X[tg == "non-targeting"].mean(0)
    D = np.vstack([X[tg == p].mean(0) - cm for p in perts])
    return D.astype(np.float64), perts, genes


def main(pred_h5ad: str) -> None:
    t0 = time.monotonic()
    r, c, emit_c = E.load_reference()
    D, perts, genes = state_effect_matrix(pred_h5ad)
    assert list(perts) == list(r["targets"]), "H1 target mismatch"
    assert genes == list(r["genes"]), "gene order mismatch vs reference"
    effects = R.map_axis(D, np.array(genes), r["genes"])
    assert np.array_equal(np.array(genes), r["genes"])

    b = R.builder()
    pick_rng, emit_rng = b.make_rngs(0)
    sums, jks, ps, qs, lfcs = [], [], [], [], []
    for i, target in enumerate(r["targets"]):
        profile = b.profile_from_effect(effects[i], r["b_ctrl"], min_cpm=10.)
        picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
        spread = b.emission_spread_vector(
            effects[i], r["b_ctrl"], min_cpm=10., tail_spread=1., keep_top=8000,
            keep_top_frac=None, keep_spread=.25, rank_score=None, spread_map="block",
            spread_lo=0., spread_hi=.60)
        gate = b.emission_rank_axis(effects[i], r["b_ctrl"], min_cpm=10.)[0]
        block, _ = rc_emit.rc_hybrid_block(
            profile, r["libraries"][picks], .1, emit_rng, spread=spread,
            overdispersion=4., control=emit_c, key=np.where(gate, effects[i], 0.),
            jitter_seed=1000 + i, call_set_first=True,
            call_set=rc_emit.CallSet(400, 2200, 1.9), emit_counts=emit_counts)
        s, j, p, q, l = E.block_stats(block, c)
        sums.append(s)
        jks.append(j)
        ps.append(p)
        qs.append(q)
        lfcs.append(l)
        if (i + 1) % 10 == 0:
            print(f"state_h1 {i + 1}/85 ({time.monotonic() - t0:.1f}s)", flush=True)
    raw, rows, _ = E.score(r, c, sums, jks, ps, qs, lfcs)
    pl.DataFrame(rows).write_csv(WORK / "per-target-state_h1.csv")

    # blend cosines on the keep-universe genes, log1p-delta space
    keep = np.asarray(c.keep)
    with np.load(BASE / "work/effects.npz", allow_pickle=False) as z:
        cap = R.map_axis(np.asarray(z["cap15"]), np.array([str(g) for g in z["genes"]]), r["genes"])
    zt = np.load(WORK / "truth_h1.npz", allow_pickle=False)
    assert [str(t) for t in zt["perts"]] == list(r["targets"])
    assert [str(g) for g in zt["genes"]] == list(r["genes"])
    T = np.asarray(zt["delta"], dtype=np.float64)[:, keep]
    S = effects[:, keep]
    C = np.asarray(cap, dtype=np.float64)[:, keep]

    def cos(a, b):
        return float((a @ b) / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-12))

    sn = S / np.maximum(np.linalg.norm(S, axis=1, keepdims=True), 1e-12)
    cn = C / np.maximum(np.linalg.norm(C, axis=1, keepdims=True), 1e-12)
    blend = (sn + cn)
    blend /= np.maximum(np.linalg.norm(blend, axis=1, keepdims=True), 1e-12)
    per = [{"target": t, "cos_state_cap15": cos(s, x),
            "cos_state_truth": cos(s, y), "cos_cap15_truth": cos(x, y),
            "cos_blend50_truth": cos(b_, y)}
           for t, s, x, y, b_ in zip(r["targets"], S, C, T, blend)]
    out = {"raw": {m: float(raw[m]) for m in raw},
           "cap15_raw": json.loads((BASE / "work/metrics-cap15.json").read_text())["raw"],
           "mean": {k: float(np.mean([d[k] for d in per])) for k in per[0] if k != "target"},
           "per_target": per, "seconds": time.monotonic() - t0}
    (WORK / "metrics-state_h1.json").write_text(json.dumps(out, indent=2))
    print(json.dumps({"raw": out["raw"], "mean_cosines": out["mean"]}, indent=2))


if __name__ == "__main__":
    main(sys.argv[1])
