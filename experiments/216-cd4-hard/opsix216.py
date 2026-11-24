"""216 — price the surviving operator-refit arm on `nmae` and the ordering half.

`op216.py` prices the refit only on `pds` and `mse`.  §6.0b measured that bringing a library
direction up to full delivered amplitude is exactly what destroys `nmae`, so no library arm
may be believed on `pds`+`mse` alone.  This re-fits the ONE operator variant that survived
(`beta = natural`, the row-count weight) and runs 188's exact identity plus the top-89
ordering precision on it, against the same a15 reference `six216.py` uses.

One eigendecomposition per fold (~6 min each, ~8 GB peak); the `equal` variant is not re-fit
because `op216.py` measured it losing on both endpoints.

    uv run python opsix216.py
"""
from __future__ import annotations

import gc
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "experiments/203-reopen-125"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, key_precision, topn, umd  # noqa: E402
from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   deliver)
from form216 import build_final, delivered_norm, solve_alpha_delivered, votes_for  # noqa: E402
from op216 import RefitOperator, basis_rows_for  # noqa: E402
from sweep216 import ALPHA, GAMMA, LAM  # noqa: E402

NMAE_SPAN = 0.60
# (name, use refit operator?, library weight, residual f, solve alpha on delivered norm?)
ARMS = [
    ("a15 f=0.10 SHIPPED", False, 0.0, 0.10, False),
    ("res f=0.34 INCUMBENT", False, 0.0, 0.34, False),
    ("refit-nat res f=0.34", True, 0.0, 0.34, False),
    ("refit-nat lib1 @adel", True, 1.0, 0.0, True),
    ("refit-nat lib0.5 @adel", True, 0.5, 0.0, True),
]


def main() -> None:
    t0 = time.time()
    from scipy.linalg import eigh as sp_eigh

    acc: dict[tuple[str, str], list] = {}
    prec: dict[str, list] = {}
    alphas: dict[str, list] = {}
    norms: dict[str, list] = {}

    for held in FOLDS:
        fold = build_fold(held)
        truth_h = np.asarray(fold["truth"], dtype=np.float64)
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        cd4e = np.asarray(cd4.mean_effect, dtype=np.float64)
        carry = np.asarray(cd4.target_condition_count) == 3
        votes = votes_for(held, fold["targets"])
        wcols = np.asarray(fold["wcols"])
        n_basis = basis_rows_for(held)

        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)

        bA = plolo.challenge_bctrl("A")
        eA = np.expm1(bA)
        gA = (1e6 * eA / float(eA.sum()) > plolo.MIN_CPM) & (bA > 0)
        sc = fold["scored"]
        cpm = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
        LR = to_lfc(truth_h, fold["b_ctrl"], gA)[:, sc]
        pools = {"wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm, 0))[None, :], POOL_N),
                 "lfc89": topn(np.abs(LR), POOL_N),
                 "h89": topn(np.abs(truth_h[:, sc]), POOL_N)}
        Ts = truth_h[:, sc]

        inc_final, _ = build_final(pre, ALPHA, op, fold, cd4, 0.10)
        target_norm = delivered_norm(inc_final, fold)
        del inc_final
        gc.collect()

        def run_arm(name, o, lw, rf, solve):
            sel = carry & (votes > 0)
            pre4 = pre.copy()
            if lw > 0:
                pre4[sel] = ((pre[sel] * votes[sel, None] + lw * cd4e[sel])
                             / (votes[sel, None] + lw))
            a = ALPHA
            if solve:
                a, _ = solve_alpha_delivered(pre4, o, fold, cd4, rf, target_norm)
            fin, _ = build_final(pre4, a, o, fold, cd4, rf)
            H = deliver(fin, bA, gA, None)
            LP = to_lfc(H, bA, gA)[:, sc]
            for pk, mm in pools.items():
                acc.setdefault((name, pk), []).append(umd(LP, LR, mm))
            prec.setdefault(name, []).append(
                key_precision(np.abs(H[:, sc]), H[:, sc], Ts, n=POOL_N))
            alphas.setdefault(name, []).append(a)
            norms.setdefault(name, []).append(
                float(np.median(np.linalg.norm(H[:, wcols], axis=1))))
            print(f"  [{held}/{name}] a {a:.4f}", flush=True)
            del fin, H, LP, pre4
            gc.collect()

        # ORDER MATTERS: the refit destroys `op.evecs` IN PLACE (the spectral
        # reconstruction scales the eigenvector matrix by sqrt(ev) without a copy, which is
        # what keeps the peak under 8 GB).  A first version of this script built the refit
        # first and then scored the two SHIPPED-operator reference arms through the wrecked
        # operator; its a15 row read `nmae` 1.05241 against `six216.py`'s 0.98476, which is
        # how the bug was caught.  Score every shipped-operator arm BEFORE touching `op`.
        for name, use_refit, lw, rf, solve in ARMS:
            if not use_refit:
                run_arm(name, op, lw, rf, solve)

        # --------------------------------------------- the refit-natural operator
        ev = np.asarray(op.evals, dtype=np.float64)
        V = np.asarray(op.evecs, dtype=np.float64)
        del op
        gc.collect()
        V *= np.sqrt(np.maximum(ev, 0.0))
        C = V @ V.T
        del V
        gc.collect()
        trC = float(np.trace(C))
        Bc = np.ascontiguousarray(
            np.where(np.asarray(cd4.gene_votes)[carry] > 0, cd4e[carry], 0.0)[:, wcols])
        m = int(Bc.shape[0])
        trC4 = float((Bc * Bc).sum() / m)
        beta = m / (n_basis + m)
        Bs = Bc * np.sqrt(beta / (m * trC4))
        Cp = np.array(C, order="F")
        Cp *= (1.0 - beta) / trC
        Cp += Bs.T @ Bs
        Cp *= trC
        del C, Bs, Bc, ev
        gc.collect()
        evp, Vp = sp_eigh(Cp, overwrite_a=True, driver="evd", check_finite=False)
        del Cp
        gc.collect()
        evp = np.maximum(evp, evp.max() * 1e-12)
        newop = RefitOperator(evp, Vp)
        print(f"[{held}] refit-natural beta {beta:.5f} ready ({time.time()-t0:.0f}s)",
              flush=True)

        for name, use_refit, lw, rf, solve in ARMS:
            if use_refit:
                run_arm(name, newop, lw, rf, solve)
        del newop, evp, Vp, pre, cd4, truth_h, fold, LR, Ts, pools
        gc.collect()

    base = ARMS[0][0]
    inc = ARMS[1][0]
    ref = {pk: np.concatenate([x[2] for x in acc[(base, pk)]])
           for pk in ("wald89", "lfc89", "h89")}
    b0 = {}
    for pk in ref:
        U = np.concatenate([x[0] for x in acc[(base, pk)]])
        M = np.concatenate([x[1] for x in acc[(base, pk)]])
        b0[pk] = float(np.mean(1 + (U - 2 * M) / ref[pk]))
    p0 = float(np.mean(prec[base]))
    a15n = float(np.mean(norms[base]))

    out = {}
    print(f"\n{'arm':<26}{'alpha':>8}{'dlv/a15':>9}{'nmae wald89':>13}{'d vs a15':>11}"
          f"{'ovr nmae':>11}{'top89 prec':>12}{'d prec':>10}{'nmae lfc89':>12}"
          f"{'nmae h89':>11}")
    for name, *_ in ARMS:
        vals = {}
        for pk in ref:
            U = np.concatenate([x[0] for x in acc[(name, pk)]])
            M = np.concatenate([x[1] for x in acc[(name, pk)]])
            vals[pk] = float(np.mean(1 + (U - 2 * M) / ref[pk]))
        nm = vals["wald89"]
        out[name] = dict(nmae_wald89=nm, d_nmae_vs_a15=nm - b0["wald89"],
                         overall_nmae_vs_a15=-(nm - b0["wald89"]) / NMAE_SPAN / 6.0,
                         top89_precision=float(np.mean(prec[name])),
                         d_precision_vs_a15=float(np.mean(prec[name])) - p0,
                         alpha=float(np.mean(alphas[name])),
                         delivered_norm_ratio=float(np.mean(norms[name])) / a15n,
                         nmae_lfc89=vals["lfc89"], nmae_h89=vals["h89"])
        o = out[name]
        print(f"{name:<26}{o['alpha']:>8.4f}{o['delivered_norm_ratio']:>9.4f}"
              f"{nm:>13.5f}{o['d_nmae_vs_a15']:>+11.5f}{o['overall_nmae_vs_a15']:>+11.6f}"
              f"{o['top89_precision']:>12.5f}{o['d_precision_vs_a15']:>+10.5f}"
              f"{o['nmae_lfc89']:>12.5f}{o['nmae_h89']:>11.5f}")

    ni = out[inc]["nmae_wald89"]
    print(f"\nagainst the incumbent `{inc}`:")
    for name, *_ in ARMS:
        o = out[name]
        o["d_nmae_vs_incumbent"] = o["nmae_wald89"] - ni
        o["overall_nmae_vs_incumbent"] = -(o["nmae_wald89"] - ni) / NMAE_SPAN / 6.0
        print(f"  {name:<26} d nmae {o['d_nmae_vs_incumbent']:>+9.5f}   "
              f"ovr {o['overall_nmae_vs_incumbent']:>+9.6f}   "
              f"d prec {o['top89_precision']-out[inc]['top89_precision']:>+8.5f}")
    (WORK / "opsix216.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {WORK/'opsix216.json'} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
