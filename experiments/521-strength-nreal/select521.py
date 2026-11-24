"""521 — selection-aware test of the rank-policy grid: the best cell of RANK(q, K) (8 q x 12 K) over the four strength
features is compared with the SAME maximum under 200 joint permutations of the target labels (every feature permuted by
one permutation, so their correlation is kept). A per-cell bootstrap CI conditional on the selected cell ignores the
selection; this null does not.
    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/select521.py
"""
from __future__ import annotations

import json
import sys

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

import lib521 as L  # noqa: E402
from policy521 import KS, QS, Table, rank_policy, value  # noqa: E402

NPERM = 200


def grid_max(T, reg, var, feats):
    R0 = T.rows[f"{reg}|base"]
    best, arg = -np.inf, None
    for fn, z in feats.items():
        for q in QS:
            for K in KS:
                v = value(T.assemble(reg, "headcpm", var, rank_policy(z, q, K)), R0)
                if v > best:
                    best, arg = v, (fn, q, K)
    return best, arg


def main():
    T = Table(0)
    P = np.load(L.WORK / "pred521_h1.npz")
    F = {"S4": P["S4__S3_p_vip_thr4"], "S3": P["S3__k562phct116phek293t_nreal521_thr4"], "S4avg": P["S4_thr-avg"],
         "K562": np.where(np.isfinite(P["k562_t4"]), P["k562_t4"], np.nanmedian(P["k562_t4"]))}
    out = {}
    for reg in ("1.0", "0.5", "0.35", "cmp"):
        for var in ("b", "a"):
            real, arg = grid_max(T, reg, var, F)
            rng = np.random.default_rng(52100)
            nul = []
            for _ in range(NPERM):
                perm = rng.permutation(85)
                nul.append(grid_max(T, reg, var, {k: v[perm] for k, v in F.items()})[0])
            nul = np.asarray(nul)
            p = float(np.mean(nul >= real))
            out[f"{reg}|{var}"] = dict(real_max=real, arg=list(map(str, arg)), null_mean=float(nul.mean()),
                                      null_q95=float(np.quantile(nul, 0.95)), p=p)
            print(f"regime {reg:4s} variant {var}: best cell {arg} {real:+.5f}; permuted-label max mean {nul.mean():+.5f}, "
                  f"95th pct {np.quantile(nul, 0.95):+.5f}; P(null max >= real) {p:.3f}", flush=True)
    (L.WORK / "select521.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
