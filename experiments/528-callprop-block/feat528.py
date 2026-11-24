#!/usr/bin/env python3
"""528 — per-gene control features of each Flex context (H1 = h1-test's 18,400 controls; A/B/C = the challenge controls),
the inputs of a call-propensity block key. Memory-bounded (517's CSRStream, fully preallocated).

Per context, control cells only:
  pooled_cpm  = 1e6 * sum_i c_ig / sum_i L_i           (517's `mean_cpm`)
  mean_cpm    = mean_i 1e6 * c_ig / L_i                (the scorer's / emitter's `ref_lin` on the 5-CPM universe)
  det         = share of control cells with c_ig > 0
  exp_det     = 1 - mean_i exp(-p_g L_i), 10 library bins (517's `excess_dropout` rule)
  x           = exp_det - det                         (excess dropout; > 0 = fewer detections than Poisson sampling)
  keep        = mean_cpm > 5                          (the scorer's DE universe)
Known answer: H1's x / det / exp_det / pooled_cpm equal 517's `excess_dropout(h1-test)` exactly, and H1's `keep` and
`mean_cpm[keep]` equal 323's reference.npz `keep` / `emit_ref_lin` (the emitter's ranking key) to float precision.

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/feat528.py [H1 A B C]
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
import h5py  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/517-flex-programme"))
from stream517 import CSRStream  # noqa: E402  (read-only import)

PATHS = {"H1": ROOT / "data/harness/h1-test.h5ad", **{c: ROOT / f"data/challenge/context_{c}.h5ad" for c in "ABC"}}


def cat(n):
    """anndata string column: categorical, nullable-string array or plain."""
    if isinstance(n, h5py.Group):
        if "categories" in n:
            return n["categories"][:].astype(str)[n["codes"][:]]
        return n["values"][:].astype(str)
    return n[:].astype(str)


def gene_names(f):
    var = f["var"]
    key = var.attrs.get("_index", "_index")
    return cat(var[key])


def features(path, nb=10):
    it = CSRStream(path)
    tv = cat(it.f["obs/target_gene"])
    genes = gene_names(it.f)
    rows = np.flatnonzero(tv == "non-targeting")
    lib = it.libraries()
    lab = np.full(it.n, -1, np.int64)
    lab[rows] = 0
    counts = it.group_sums(lab, 1)[0]
    dets = it.group_sums(lab, 1, binary=True)[0]
    cpms = it.group_sums(lab, 1, cpm=True, lib=lib)[0]
    it.close()
    n = rows.size
    p = counts / counts.sum()
    L = lib[rows]
    edges = np.quantile(L, np.linspace(0, 1, nb + 1))
    b = np.clip(np.searchsorted(edges, L, side="right") - 1, 0, nb - 1)
    Lb = np.array([L[b == k].mean() for k in range(nb)])
    nbk = np.bincount(b, minlength=nb) / n
    exp_det = 1.0 - (nbk[None, :] * np.exp(-p[:, None] * Lb[None, :])).sum(1)
    det = dets / n
    mean_cpm = cpms / n
    return dict(genes=genes, pooled_cpm=1e6 * p, mean_cpm=mean_cpm, det=det, exp_det=exp_det, x=exp_det - det,
                keep=mean_cpm > 5.0, n_ctrl=np.int64(n), lib_median=np.float64(np.median(L)))


def main(ctxs):
    out = {}
    for c in ctxs:
        t0 = time.time()
        F = features(PATHS[c])
        np.savez(HERE / f"work/feat528_{c}.npz", **F)
        e = F["keep"]
        out[c] = dict(n_ctrl=int(F["n_ctrl"]), n_keep=int(e.sum()), lib_median=float(F["lib_median"]),
                      x_median_keep=float(np.median(F["x"][e])), x_q90_keep=float(np.quantile(F["x"][e], .9)),
                      frac_x_gt_0_05=float((F["x"][e] > 0.05).mean()), seconds=time.time() - t0)
        print(c, json.dumps(out[c]), flush=True)
        if c == "H1":
            # KA 1: 517's excess_dropout on h1-test (its own function, read-only import)
            from dropout517 import excess_dropout  # noqa: PLC0415
            H = excess_dropout(PATHS["H1"])
            ka = {k: float(np.max(np.abs(H[k] - F[{"mean_cpm": "pooled_cpm"}.get(k, k)]))) for k in ("x", "det", "exp_det", "mean_cpm")}
            # KA 2: the reference's gate and the emitter's ranking key
            z = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
            rg = z["genes"].astype(str)
            assert np.array_equal(rg, F["genes"]), "h1-test gene axis differs from 323's reference axis"
            keep_eq = bool(np.array_equal(np.asarray(z["keep"], bool), F["keep"]))
            d_lin = float(np.max(np.abs(np.asarray(z["emit_ref_lin"], float) - F["mean_cpm"][F["keep"]]) /
                                 np.asarray(z["emit_ref_lin"], float))) if keep_eq else float("nan")
            rank_eq = bool(np.array_equal(np.argsort(-np.asarray(z["emit_ref_lin"], float), kind="stable"),
                                          np.argsort(-F["mean_cpm"][F["keep"]], kind="stable"))) if keep_eq else False
            out["KA_H1"] = dict(vs_517_maxabs=ka, keep_equal=keep_eq, emit_ref_lin_max_rel=d_lin, emit_rank_equal=rank_eq)
            print("KA H1:", json.dumps(out["KA_H1"]), flush=True)
    prev = json.loads((HERE / "work/feat528.json").read_text()) if (HERE / "work/feat528.json").exists() else {}
    prev.update(out)
    (HERE / "work/feat528.json").write_text(json.dumps(prev, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:] or ["H1", "A", "B", "C"])
