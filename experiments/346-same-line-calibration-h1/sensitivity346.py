#!/usr/bin/env python3
"""346 step 7 — how many of the recipient's own training targets the sign gain needs.

Refits the same-line full ridge on random subsets of the 1,965 Pan H1 training targets and
reports the design-stage |lfc|-weighted sign accuracy of the beta=1 integration for each size.
This is the number the PI needs to decide whether context A's non-panel bank is large enough.
Design stage only: no emission, no scorer call. Writes only into this experiment's work/.
"""
from __future__ import annotations

import gc
import json
import os
import resource
import sys
import time
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/334-text-conditioned-transfer"))
sys.path.insert(0, str(ROOT / "experiments/341-known-line-calibration"))
import core334 as P      # noqa: E402
import core341 as C      # noqa: E402
import fit346 as F       # noqa: E402
import effects346 as X   # noqa: E402

SIZES = (125, 250, 500, 1000, 1965)
SEEDS = (0, 1, 2)


def main() -> None:
    clock = time.monotonic()
    axes = F.build_axes()
    direction, basal, rms, _ = F.assemble(axes, poison=False)
    n_train = len(axes["training"])
    test = np.arange(n_train, direction.shape[1])
    donors = [0, 1, 2]
    query = direction[donors].mean(0)[test]
    predictions = {}
    for size in SIZES:
        for seed in SEEDS:
            if size == n_train and seed != SEEDS[0]:
                continue           # the full set has no sampling variation
            rng = np.random.default_rng(3460910 + 1000 * seed + size)
            sub = np.sort(rng.choice(n_train, size=min(size, n_train), replace=False))
            inner = P.source_blocks(direction, basal, rms, 3, sub)
            x = np.concatenate([f[..., 0] for _, f, _ in inner])
            y = np.tile(direction[3, sub], (len(donors), 1))
            model = C.fit(x, y)
            predictions[f"n{size}_s{seed}"] = C.predict(model, query)
            del inner, x, y, model
            gc.collect()
            print(f"fitted n={size} seed={seed} ({time.monotonic()-clock:.0f}s)", flush=True)
    del direction, basal, rms
    gc.collect()

    with np.load(BASE / "work/effects.npz", allow_pickle=False) as z:
        genes = [str(g) for g in z["genes"]]
        targets = [str(t) for t in z["targets"]]
        cap15 = np.asarray(z["cap15"], dtype=np.float64)
    gi = {g: i for i, g in enumerate(genes)}
    cols = np.array([gi[g] for g in axes["covered"]], dtype=np.int64)
    ti = {t: i for i, t in enumerate(targets)}
    rows = np.array([ti[t] for t in axes["test"]], dtype=np.int64)
    arms = {"cap15": cap15}
    for key, value in predictions.items():
        new = cap15.copy()
        for j, i in enumerate(rows):
            new[i, cols] = X.blend(cap15[i, cols], value[j], 1.0)
        arms[key] = new

    with np.load(BASE / "work/reference.npz", allow_pickle=False, mmap_mode="r") as z:
        genes_all = np.asarray(z["genes"]).astype(str)
        keep = np.asarray(z["keep"], dtype=bool)
        ref_targets = np.asarray(z["targets"]).astype(str)
    dgenes = genes_all[keep]
    dgi = {g: i for i, g in enumerate(dgenes)}
    eg = {g: i for i, g in enumerate(genes)}
    col = np.array([eg.get(g, -1) for g in dgenes])
    miss = col < 0
    mapped = {}
    for a, e in arms.items():
        m = e[:, np.where(miss, 0, col)]
        m[:, miss] = 0.
        mapped[a] = m
    del arms
    gc.collect()
    ref = pl.read_parquet(BASE / "work/reference_de.parquet")
    values = {a: [] for a in mapped}
    for tix, t in enumerate(ref_targets):
        sub = ref.filter(pl.col("target") == t)
        g = sub["feature"].to_numpy().astype(str)
        q = sub["p_adj"].to_numpy()
        l = sub["log2_fold_change"].to_numpy().astype(float)
        gate = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
        if gate.sum() < 10:
            continue
        idx = np.array([dgi[x] for x in g[gate]])
        l = l[gate]
        w = np.abs(l)
        D = w.sum()
        for a in mapped:
            se = np.sign(mapped[a][tix, idx]) == np.sign(l)
            values[a].append(float((w * se).sum() / D))
    values = {a: np.asarray(v) for a, v in values.items()}
    out = {"n_eligible_targets": int(len(values["cap15"])),
           "cap15_design_weighted_sign": float(values["cap15"].mean()),
           "by_training_size": {}}
    for size in SIZES:
        gains = [float((values[k] - values["cap15"]).mean())
                 for k in values if k.startswith(f"n{size}_")]
        out["by_training_size"][str(size)] = {
            "n_fits": len(gains), "mean_gain": float(np.mean(gains)),
            "gains": gains,
            "mean_design_weighted_sign": float(np.mean(
                [float(values[k].mean()) for k in values if k.startswith(f"n{size}_")])),
        }
    out["seconds"] = time.monotonic() - clock
    out["max_rss_bytes_macos"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    out["scope"] = ("design-stage |lfc|-weighted sign accuracy after the beta=1 integration; "
                    "not a scored member")
    (WORK / "sensitivity.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
