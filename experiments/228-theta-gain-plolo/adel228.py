"""228 Phase 4 — THE CRUX: the gain at AMPLITUDE-MATCHED alpha.

Phase 2 preserved the median row norm of the EFFECT matrix over `wcols`, which is
production's own convention and makes the arm amplitude-neutral in the space `VCC_ALPHA`
acts in.  It is NOT amplitude-neutral in the space the SCORER reads: `profile_from_effect`
is nonlinear, and Phase 3 measured the DELIVERED median row norm up **+7.4%** for
`gain_reg@0.25` on the three cancer folds -- with `nmae` degrading in near-lockstep across
every arm and every null.

211's standing rule says `claim` is at 100.0000% of its cap, so that +7.4% is amplitude we do
not have.  This script therefore re-solves `VCC_ALPHA` per arm by 216's `@adel` policy --
a secant solve on `log(delivered median row norm)` against the SHIPPED STACK's own value --
and re-reads all six members there.  The difference between the two readings is exactly
"how much of the gain is amplitude", which is what 224 §10 (ii) asked for.

    uv run python adel228.py --arms gain_reg@0.25 ...
"""
from __future__ import annotations

import argparse
import hashlib
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
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, key_precision, topn, umd  # noqa: E402
from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   deliver, solve_fraction)

sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
from resid218 import ALPHA, GAMMA, LAM, score_effects, source_residual  # noqa: E402
from compose218 import merged, solve_fraction_carrying, union_rows  # noqa: E402

from arms228 import CD4_F, UNION, UNION_G, apply_gain, multiplier  # noqa: E402
from ctrl228 import control_block  # noqa: E402
from six228 import (CELLS, MEDIAN_LIBRARY, N_REALS, block_digest,  # noqa: E402
                    build_final, call_set)

TOL, MAXIT = 2e-4, 10


def delivered_norm(final, fold, bctx, gate) -> float:
    vals = []
    for c in ("A", "B", "C"):
        dl = deliver(final, bctx[c], gate[c], "minnorm")
        vals.append(float(np.median(np.linalg.norm(dl[:, fold["scored"]], axis=1))))
    return float(np.mean(vals))


def solve_alpha(make, fold, bctx, gate, target, a0=1.0):
    xs, ys = [a0, a0 * 0.9], []
    for a in xs:
        f = make(a)
        ys.append(np.log(delivered_norm(f, fold, bctx, gate)) - np.log(target))
        del f
    for _ in range(MAXIT):
        if abs(ys[-1]) < TOL:
            break
        den = ys[-1] - ys[-2]
        if den == 0:
            break
        a = float(np.clip(xs[-1] - ys[-1] * (xs[-1] - xs[-2]) / den, 0.2, 5.0))
        f = make(a)
        xs.append(a)
        ys.append(np.log(delivered_norm(f, fold, bctx, gate)) - np.log(target))
        del f
    return float(xs[-1]), float(np.exp(ys[-1]))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="*", default=["gain_reg@0.25"])
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    args = ap.parse_args()
    out_path = WORK / "adel228.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    ctrl = {c: control_block(c) for c in ("A", "B", "C")}
    bctx = {c: plolo.challenge_bctrl(c) for c in ("A", "B", "C")}
    gate = {}
    for c in ("A", "B", "C"):
        e = np.expm1(bctx[c])
        gate[c] = (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (bctx[c] > 0)

    for held in args.folds:
        t0 = time.time()
        fold = build_fold(held)
        wcols = np.asarray(fold["wcols"])
        sc = np.asarray(fold["scored"])
        truth_h = np.asarray(fold["truth"], dtype=np.float64)
        Ts = truth_h[:, sc]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        h1 = union_rows([source_residual(n, fold["targets"], "gene_norm") for n in UNION])
        op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                      if FOLDS[held]["op_kind"] == "production" else Path(FOLDS[held]["op"]))
        z = np.load(WORK / f"gain-{held}.npz", allow_pickle=True)
        live = np.asarray(z["live"], dtype=bool)

        cpm_sc = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
        LR = to_lfc(truth_h, fold["b_ctrl"], gate["A"])[:, sc]
        pools = {"wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm_sc, 0)), POOL_N),
                 "lfc89": topn(np.abs(LR), POOL_N),
                 "h89": topn(np.abs(Ts), POOL_N)}
        real_de = {n: topn(np.abs(Ts), n) for n in N_REALS}

        base_final = build_final(fold, op, pre, cd4, h1, wcols, None, 1.0)
        target = delivered_norm(base_final, fold, bctx, gate)
        del base_final
        print(f"[{held}] shipped-stack delivered median row norm {target:.6f}", flush=True)

        def measure(final, arm, a, resid):
            H_A = deliver(final, bctx["A"], gate["A"], None)
            LP = to_lfc(H_A, bctx["A"], gate["A"])[:, sc]
            nm = {k: umd(LP, LR, p) for k, p in pools.items()}
            prec89 = key_precision(np.abs(H_A[:, sc]), H_A[:, sc], Ts, n=POOL_N)
            fid, jac, digests = [], {n: [] for n in N_REALS}, []
            for c in ("A", "B", "C"):
                keep, ref_lin = ctrl[c]
                ok = tot = 0
                calls = np.zeros((final.shape[0], len(sc)), dtype=bool)
                for i in range(final.shape[0]):
                    p = plolo.profile_from_effect(final[i], bctx[c], gate[c], "minnorm")
                    digests.append(block_digest(p))
                    h_row = np.log1p(plolo.TS * p) - bctx[c]
                    call, _ = call_set(final[i], gate[c], keep, ref_lin)
                    cs = call[sc]
                    calls[i] = cs
                    sp, st = np.sign(h_row[sc]), np.sign(Ts[i])
                    li = cs & (st != 0) & (sp != 0)
                    ok += int((sp[li] == st[li]).sum()); tot += int(li.sum())
                fid.append(ok / max(tot, 1))
                for n in N_REALS:
                    calls_n = (calls & real_de[n]).sum(1)
                    un = (calls | real_de[n]).sum(1)
                    jac[n].append(float(np.mean(calls_n / np.maximum(un, 1))))
            e = score_effects(final, fold, truth_h[:, sc])
            e = {k: v for k, v in e.items() if not k.startswith("_")}
            return dict(
                fold=held, arm=arm, a=float(a), alpha_effective=float(ALPHA * a),
                delivered_norm=float(delivered_norm(final, fold, bctx, gate)),
                delivered_norm_ratio=float(
                    delivered_norm(final, fold, bctx, gate) / target),
                resid=float(resid),
                nmae_pools={k: float(np.mean(1 + (u - 2 * m) / d))
                            for k, (u, m, d) in nm.items()},
                nmae_U={k: float(np.mean(v[0] / v[2])) for k, v in nm.items()},
                nmae_M={k: float(np.mean(v[1] / v[2])) for k, v in nm.items()},
                reach_prec89=float(prec89), fid_precision=float(np.mean(fid)),
                jac={str(n): float(np.mean(v)) for n, v in jac.items()},
                digest=hashlib.sha256("".join(digests).encode()).hexdigest()[:16],
                **e)

        for arm in ["BASE"] + list(args.arms):
            if (held, arm) in done:
                continue
            if arm == "BASE":
                m = None
            else:
                name, _, tt = arm.partition("@")
                m = multiplier(np.asarray(z[name], dtype=np.float64), live, float(tt))

            def make(a, m=m):
                return build_final(fold, op, pre, cd4, h1, wcols, m, a)

            a, resid = (1.0, 1.0) if arm == "BASE" else solve_alpha(
                make, fold, bctx, gate, target)
            final = make(a)
            row = measure(final, arm, a, resid)
            rows.append(row)
            out_path.write_text(json.dumps(rows, indent=2))
            print(f"[{held}/{arm}] a* {a:.5f} (alpha {ALPHA * a:.5f})  "
                  f"norm ratio {row['delivered_norm_ratio']:.5f}  "
                  f"cos_agg {row['cos_agg']:.6f}  pds {row['pds']:.6f}  "
                  f"nmae {row['nmae_pools']['wald89']:.5f}  prec89 {row['reach_prec89']:.5f}"
                  f"  ({time.time() - t0:.0f}s)", flush=True)
            del final
        del op

    out_path.write_text(json.dumps(rows, indent=2))
    print(f"\n-> {out_path}")


if __name__ == "__main__":
    main()
