"""318 step 1c — WHY the parity arm is not shape-only as delivered.

316's KA-3 certified the arm "exactly shape-only": max |delivered row norm / incumbent - 1| =
4.44e-16. That statement is made on the PRE-WHITEN effect matrix, one stage before the
operator. `profile318.py`/`family318.py` then measured what the scorer actually reads and
found the parity arm's per-target delivered magnitude ratio at sd(log) = 0.0485, range
0.923-1.307 -- LARGER than the shipped-vs-parity difference the mandate calls "the
distortion" (0.0329).

This script isolates the stage responsible, on the production operator, with no donor model:

    x_t          = 170's pre-whiten challenge candidate rows (`pre.npz: cand_chal`)
    y_t          = x_t rotated by a fixed angle into a random orthogonal direction and then
                   RENORMALISED to ||x_t||   -- i.e. exactly what a norm match delivers:
                   a different direction at an identical pre-whiten norm
    g(t)         = ||W x_t|| / ||x_t||       -- the operator's per-row gain
    delivered(t) = ||W y_t|| / ||W x_t||     -- what a norm-matched rotation delivers

If the pre-whiten norm match controlled the delivered magnitude, `delivered(t)` would be 1 for
every target. The spread of `delivered` at the parity arm's own measured rotation is the
prediction for the delivered redistribution 316's KA-3 cannot see.

Peak RSS is dominated by one 17,942 x 17,942 float64 eigenbasis (2.6 GiB) plus the GEMM.
"""
from __future__ import annotations

import importlib.util
import json
import os
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")
os.environ.setdefault("VCC_WHITEN_GAMMA", "0.65")
os.environ.setdefault("VCC_WHITEN_LAM", "0.0")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
sys.path.insert(0, str(REPO / "src"))


def rss() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 ** 3) if sys.platform == "darwin" else r / (1024 ** 2)


def main() -> None:
    t0 = time.time()
    spec = importlib.util.spec_from_file_location(
        "prod_builder_318", REPO / "scripts/make_submission_pm_cd4.py")
    B = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_318"] = B
    spec.loader.exec_module(B)
    B.assert_frozen_contract()

    pre = np.load(REPO / "experiments/170-gamma-mse/work/pre.npz", allow_pickle=True)
    cols = pre["operator_columns"].astype(np.int64)
    X = np.asarray(pre["cand_chal"], dtype=np.float64)
    whit = B.Whitener.load(REPO / "data/donors/whitener_union_pruned.npz")
    print(f"[318] X {X.shape}, operator {whit.evecs.shape}, RSS {rss():.2f} GiB", flush=True)

    # the measured rotations: parity and shipped against the incumbent, delivered
    fam = json.loads((HERE / "work/predict318.json").read_text())
    cosp = float(np.median(fam["cos_parity"]))
    coss = float(np.median(fam["cos_shipped"]))

    WX, _ = B.apply_v9_whitening(X, whit, cols, gamma=0.65, lam=0.0)
    nx = np.linalg.norm(X[:, cols], axis=1)
    nwx = np.linalg.norm(WX, axis=1)
    g = nwx / np.maximum(nx, 1e-300)
    print(f"[318] operator per-row gain g(t): median {np.median(g):.6f} range "
          f"[{g.min():.6f}, {g.max():.6f}] sd(log) {np.std(np.log(g)):.6f}", flush=True)

    out = {"n_targets": int(X.shape[0]),
           "gain": {"median": float(np.median(g)), "min": float(g.min()),
                    "max": float(g.max()), "sd_log": float(np.std(np.log(g)))},
           "measured_median_cos": {"parity_vs_incumbent": cosp,
                                   "shipped_vs_incumbent": coss},
           "rotations": {}}

    rng = np.random.default_rng(318)
    Xc = X[:, cols]
    for label, c in (("KA_cos_1.0_identity", 1.0),
                     ("parity_median_cos", cosp), ("shipped_median_cos", coss),
                     ("cos_0.999", 0.999), ("cos_0.99", 0.99), ("cos_0.95", 0.95)):
        tan = float(np.sqrt(max(1.0 / c ** 2 - 1.0, 0.0)))
        V = rng.standard_normal(Xc.shape)
        V -= (V * Xc).sum(1, keepdims=True) / np.maximum(
            (Xc * Xc).sum(1, keepdims=True), 1e-300) * Xc
        V *= (np.linalg.norm(Xc, axis=1) / np.maximum(
            np.linalg.norm(V, axis=1), 1e-300))[:, None]
        # a COPY of X: the panel columns outside the operator must not move, or the
        # preserve scale changes and the KA row at cos = 1 stops reading exactly 1.
        Y = X.copy()
        Y[:, cols] = Xc + tan * V
        # the norm match: identical PRE-WHITEN norm, different direction
        Y[:, cols] *= (nx / np.maximum(np.linalg.norm(Y[:, cols], axis=1), 1e-300))[:, None]
        chk = float(np.abs(np.linalg.norm(Y[:, cols], axis=1) / nx - 1.0).max())
        WY, _ = B.apply_v9_whitening(Y, whit, cols, gamma=0.65, lam=0.0)
        r = np.linalg.norm(WY, axis=1) / np.maximum(nwx, 1e-300)
        cosreal = float(np.median((Xc * Y[:, cols]).sum(1)
                                  / (nx * np.linalg.norm(Y[:, cols], axis=1))))
        out["rotations"][label] = {
            "asked_cos": c, "realised_median_cos": cosreal,
            "prewhiten_norm_match_max_dev": chk,
            "delivered_ratio": {"median": float(np.median(r)), "min": float(r.min()),
                                "max": float(r.max()), "sd_log": float(np.std(np.log(r))),
                                "whole": float(np.sqrt((np.linalg.norm(WY, axis=1) ** 2).sum()
                                                       / (nwx ** 2).sum()))}}
        print(f"[318] rotation {label:20s} cos {cosreal:.6f}: pre-whiten norms matched to "
              f"{chk:.2e}; DELIVERED norm ratio median {np.median(r):.6f} range "
              f"[{r.min():.6f}, {r.max():.6f}] sd(log) {np.std(np.log(r)):.6f}", flush=True)
        del V, Y, WY

    out["peak_rss_gib"] = rss()
    out["seconds"] = round(time.time() - t0, 1)
    (HERE / "work/mechanism318.json").write_text(json.dumps(out, indent=1))
    print(f"[318] peak RSS {rss():.2f} GiB, {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
