"""690b — a depth-matched SAME-PANEL mock D/E/F bundle ("same20k") for the AUTO_R0 end-to-end test (2026-10-15).

Arc releases D/E/F at A/B/C's depth (~20k median UMI, the 10-01 post), so the real day passes 641's depth row with no
override. 690's bundles carry Orion controls at x0.80-0.84 of that depth and need DEF_DAY_ALLOW_DEPTH=1. This bundle
matches A's library-size distribution instead, so the day's whole automated path (691's watcher with AUTO_R0=1) runs with
NO override:
  D = HCT116 (Orion NTC): a pool of 55,200 NTC cells, importance-subsampled without replacement to 18,400 cells whose
      log-library histogram follows context A's (40 bins; weight = p_A / p_pool per bin)
  E = H1 (Arc 2025 Flex controls): 690's context_E.h5ad, already at 20,186 median (copied, relabelled nothing)
  F = HEK293T (Orion NTC): as D
Panel: A/B/C's 300 targets (the SAME-PANEL branch). Selecting cells by library biases composition slightly toward large
cells; for a rehearsal of the day's mechanics that is acceptable and stated.

Usage: uv run python experiments/690-def-readiness-1015/build690b.py <out bundle dir>
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import sys
import time
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("b595", ROOT / "experiments/595-foreign-control-rehearsal/build595.py")
b595 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b595)
POOL, N, SEED = 55_200, 18_400, 690
E_SRC = ROOT / "experiments/690-def-readiness-1015/work/bundles/new/context_E.h5ad"


def lib_sizes_csr(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as h:
        X = h["X"]
        ptr = X["indptr"][:].astype(np.int64)
        data = X["data"]
        out = np.empty(len(ptr) - 1)
        step = 2000
        for a in range(0, len(ptr) - 1, step):
            b = min(len(ptr) - 1, a + step)
            blk = data[ptr[a]:ptr[b]].astype(np.float64)
            out[a:b] = np.add.reduceat(blk, ptr[a:b] - ptr[a]) if len(blk) else 0.0
        return out


def matched(m, ref_lib: np.ndarray, rng) -> tuple:
    lib = np.asarray(m.sum(1)).ravel()
    lo, hi = np.log10(max(1.0, min(lib.min(), ref_lib.min()))), np.log10(max(lib.max(), ref_lib.max())) + 1e-9
    edges = np.linspace(lo, hi, 41)
    pa = np.histogram(np.log10(ref_lib), edges)[0] / len(ref_lib)
    pp = np.histogram(np.log10(np.maximum(lib, 1)), edges)[0] / len(lib)
    b = np.clip(np.digitize(np.log10(np.maximum(lib, 1)), edges) - 1, 0, 39)
    w = np.where(pp[b] > 0, pa[b] / np.maximum(pp[b], 1e-12), 0.0)
    pick = np.sort(rng.choice(len(lib), N, replace=False, p=w / w.sum()))
    return m[pick], {"pool": int(len(lib)), "pool_lib_median": float(np.median(lib)),
                     "picked_lib_median": float(np.median(lib[pick])), "ref_lib_median": float(np.median(ref_lib))}


def main() -> None:
    out = Path(sys.argv[1]).resolve()
    if out.exists():
        raise SystemExit(f"{out} exists")
    out.mkdir(parents=True)
    rng = np.random.default_rng(SEED)
    ref = lib_sizes_csr(ROOT / "data/challenge/context_A.h5ad")
    b595.log(f"context A library: median {np.median(ref):.0f}, q25 {np.percentile(ref, 25):.0f}, q75 {np.percentile(ref, 75):.0f}")
    shutil.copy(ROOT / "data/challenge/gene_names.csv", out / "gene_names.csv")
    shutil.copy(ROOT / "data/challenge/pert_counts.csv", out / "pert_counts.csv")
    man = json.loads((ROOT / "experiments/180-def-readiness/work/def-bundle/manifest.json").read_text())
    man.update(partition="synthetic-def", panel_id="vcc2026-mock-690b-same20k")
    (out / "manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    record = {"seed": SEED, "pool": POOL, "panel": "data/challenge/pert_counts.csv", "contexts": {}}
    for ctx, line in (("D", "HCT116"), ("F", "HEK293T")):
        t0 = time.time()
        b595.N_CELLS = POOL          # read a larger pool (orion() samples N_CELLS without replacement from what it read)
        m, src = b595.orion(line)
        b595.N_CELLS = N             # write_context builds 46 x 400 ntc ids from N_CELLS
        m, info = matched(m, ref, rng)
        rec = b595.write_context(m, ctx, out / f"context_{ctx}.h5ad")
        rec.update(line=line, **src, **info, seconds=round(time.time() - t0, 1))
        record["contexts"][ctx] = rec
        b595.log(f"  {ctx} {line}: median lib {rec['lib_median']:.0f} (pool {info['pool_lib_median']:.0f}), q75/q25 {rec['q75_q25']:.3f}")
        del m
    shutil.copy(E_SRC, out / "context_E.h5ad")
    record["contexts"]["E"] = {"line": "H1", "copied_from": str(E_SRC.relative_to(ROOT)),
                               "sha256": hashlib.sha256((out / "context_E.h5ad").read_bytes()).hexdigest()}
    (out.parent / "build690b.json").write_text(json.dumps(record, indent=1, default=str) + "\n")
    b595.log(f"wrote {out}")


if __name__ == "__main__":
    main()
