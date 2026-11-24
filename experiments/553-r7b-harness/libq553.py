#!/usr/bin/env python3
"""553 — control-library quantiles of every Flex context on disk and of R7b (control cells only; per-cell totals on the
official axis = metadata-level statistics of control cells). The production emitter's classes are q25 / q50 / q75.
usage: uv run python experiments/553-r7b-harness/libq553.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/517-flex-programme"))
from stream517 import CSRStream  # noqa: E402

PATHS = {"A": ROOT / "data/challenge/context_A.h5ad", "B": ROOT / "data/challenge/context_B.h5ad",
         "C": ROOT / "data/challenge/context_C.h5ad", "R7b": ROOT / "experiments/552-r7-reserve/work/r7b/bundle/context_R.h5ad"}
out = {}
for c, p in PATHS.items():
    it = CSRStream(p)
    lib = it.libraries()
    it.close()
    q = np.quantile(lib, [0, .01, .1, .25, .5, .75, .9, 1])
    out[c] = dict(n=int(lib.size), q=[float(v) for v in q], q75_over_q25=float(q[5] / q[3]), q50_over_q25=float(q[4] / q[3]))
z = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz")
lib = z["libraries"]
q = np.quantile(lib, [0, .01, .1, .25, .5, .75, .9, 1])
out["H1"] = dict(n=int(lib.size), q=[float(v) for v in q], q75_over_q25=float(q[5] / q[3]), q50_over_q25=float(q[4] / q[3]))
(HERE / "work/libq553.json").write_text(json.dumps(out, indent=1))
print("| context | controls | min | q1 | q10 | q25 | q50 | q75 | q90 | max | q75/q25 |")
print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for c, v in out.items():
    print(f"| {c} | {v['n']} | " + " | ".join(f"{x:,.0f}" for x in v["q"]) + f" | {v['q75_over_q25']:.2f} |")
