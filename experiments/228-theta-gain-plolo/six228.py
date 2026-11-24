"""228 Phase 3 — the members PLOLO's cosine/`pds` endpoint does not carry: `nmae`, the
ordering half of `reach`, and -- the thing 224 said must not be carried at zero -- the DE
CALL SET, hence `fid` and `jac`.

THE STRUCTURAL FACT THAT MAKES `fid`/`jac` PRICEABLE HERE
--------------------------------------------------------
The shipped emitter is `rc_callset` (board entry `31oB9yWy8CSNBDEkY4iW`, +0.0084448), whose
call set is PRESCRIBED:

    head = top `RC_CALLSET_HEAD` = 400 genes by |effect|, inside the scorer's 5-CPM universe
    cpm  = top `RC_CALLSET_CPM`  = 2,200 more by the CONTROL's mean per-cell CPM
    n_pred = 2,600 exactly, enforced by `z_cap = 1.90`

So under the shipped emitter a per-gene multiplier on the effect matrix:
  * CANNOT move `n_pred` -- it is a build constant, and `min(1, n_pred/n_real)` is already
    saturated (`n_real` median 89, ours 2,600), so `fid`'s yield half is pinned at 1;
  * CANNOT move 2,200 of the 2,600 called genes -- that block is a function of the control
    cells alone;
  * can move only WHICH 400 genes are in the head, and the SIGNS the scorer reads.

That converts "`fid` and `jac` are unpriceable" (rule 9, and 216's own §3 declaration) into a
bounded, measurable quantity on this instrument.

`nmae` is 188's exact identity with the `|lfc|*sqrt(cpm)` pool at the reference's BH depth
(89 genes/row), the construction 188 validated on h1-test, with both alternative pools
carried because 188 showed the construction can invert the sign.  `reach`'s ordering half is
the |truth|-weighted sign precision of the top-89 coordinates under the production key |H|;
its yield half carries `n_real` and is unreadable on any local harness (rule 9).

DIGEST (rule 13's single-valued test).  The per-gene panel-total digest is computed on the
PLOLO frame's EXPECTED totals (400 cells x the 20,109-count median library x the delivered
profile) for every (context, target) block, and compared against the shipped stack's.

    uv run python six228.py --arms BASE gain@0.1 gain_reg@0.25 ...
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
from resid218 import ALPHA, GAMMA, LAM, source_residual  # noqa: E402
from compose218 import merged, solve_fraction_carrying, union_rows  # noqa: E402

from arms228 import CD4_F, UNION, UNION_G, apply_gain, multiplier  # noqa: E402
from ctrl228 import control_block  # noqa: E402

RC_HEAD, RC_CPM = 400, 2200
MEDIAN_LIBRARY = 20109.0        # AGENTS.md §2, the reference's median UMI/cell
CELLS = 400
NMAE_SPAN = 0.60
N_REALS = (89, 200, 400)


def call_set(effect_row, gate10, keep, ref_lin):
    """`rc_emit.CallSet(400, 2200, 1.90)` exactly: |effect| head, then control-CPM block."""
    score = np.where(gate10 & keep, np.abs(np.nan_to_num(effect_row)), -np.inf)
    head = np.argsort(-score, kind="stable")[:RC_HEAD]
    rest = np.where(keep, ref_lin, -np.inf).copy()
    rest[head] = -np.inf
    cpm = np.argsort(-rest, kind="stable")[:RC_CPM]
    call = np.zeros(len(score), dtype=bool)
    call[head] = True
    call[cpm] = True
    return call, head


def block_digest(profile: np.ndarray) -> str:
    tot = np.rint(CELLS * MEDIAN_LIBRARY * profile).astype(np.int64)
    return hashlib.sha256(tot.tobytes()).hexdigest()[:16]


def build_final(fold, op, pre, cd4, h1, wcols, m, a=1.0):
    whitened = op.apply(pre * (ALPHA * a), wcols, GAMMA, LAM)
    if m is not None:
        whitened, _ = apply_gain(whitened, m, wcols)
    frac_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
    frac_u, _ = solve_fraction_carrying(whitened, h1, UNION_G)
    final, _ = E062.add_post_whiten_residual(
        whitened, merged(cd4, frac_cd4, h1, frac_u), 1.0, wcols)
    final[:, fold["panel_positions"]] = 0.0
    return final


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="*", default=["BASE"])
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--out", default=str(WORK / "six228.json"))
    args = ap.parse_args()
    out_path = Path(args.out)
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

        for arm in args.arms:
            if (held, arm) in done:
                continue
            if arm == "BASE":
                m, a = None, 1.0
            else:
                spec, _, atxt = arm.partition("|a=")
                name, _, tt = spec.partition("@")
                m = multiplier(np.asarray(z[name], dtype=np.float64), live, float(tt))
                a = float(atxt) if atxt else 1.0
            final = build_final(fold, op, pre, cd4, h1, wcols, m, a)

            # -------- nmae (188's identity, challenge frame A, a15 render as in 207/216)
            H_A = deliver(final, bctx["A"], gate["A"], None)
            LP = to_lfc(H_A, bctx["A"], gate["A"])[:, sc]
            nm = {k: umd(LP, LR, p) for k, p in pools.items()}

            # -------- reach ordering half, production key |H|
            prec89 = key_precision(np.abs(H_A[:, sc]), H_A[:, sc], Ts, n=POOL_N)

            # -------- the call set, per challenge context, on the a23 delivered profile
            fid, jac, digests = [], {n: [] for n in N_REALS}, []
            head_ids = []
            for c in ("A", "B", "C"):
                keep, ref_lin = ctrl[c]
                ok, tot = 0, 0
                calls = np.zeros((final.shape[0], len(sc)), dtype=bool)
                hset = []
                for i in range(final.shape[0]):
                    p = plolo.profile_from_effect(final[i], bctx[c], gate[c], "minnorm")
                    digests.append(block_digest(p))
                    h_row = np.log1p(plolo.TS * p) - bctx[c]
                    call, head = call_set(final[i], gate[c], keep, ref_lin)
                    hset.append(np.sort(head))
                    cs = call[sc]
                    calls[i] = cs
                    s_pred = np.sign(h_row[sc]); s_true = np.sign(Ts[i])
                    live_i = cs & (s_true != 0) & (s_pred != 0)
                    ok += int((s_pred[live_i] == s_true[live_i]).sum())
                    tot += int(live_i.sum())
                fid.append(ok / max(tot, 1))
                head_ids.append(np.stack(hset))
                for n in N_REALS:
                    inter = (calls & real_de[n]).sum(1)
                    union = (calls | real_de[n]).sum(1)
                    jac[n].append(float(np.mean(inter / np.maximum(union, 1))))
            head_j = np.concatenate([h.ravel() for h in head_ids])

            row = dict(fold=held, arm=arm, a=a, n_targets=final.shape[0],
                       nmae_pools={k: float(np.mean(1 + (u - 2 * mm) / d))
                                   for k, (u, mm, d) in nm.items()},
                       nmae_U={k: [float(x) for x in v[0]] for k, v in nm.items()},
                       nmae_M={k: [float(x) for x in v[1]] for k, v in nm.items()},
                       nmae_D={k: [float(x) for x in v[2]] for k, v in nm.items()},
                       reach_prec89=float(prec89),
                       fid_precision=float(np.mean(fid)),
                       fid_precision_per_ctx=[float(x) for x in fid],
                       n_pred=int(RC_HEAD + RC_CPM),
                       jac={str(n): float(np.mean(v)) for n, v in jac.items()},
                       digest=hashlib.sha256("".join(digests).encode()).hexdigest()[:16],
                       block_digests=digests,
                       head_members=head_j.tolist(),
                       delivered_median_norm=float(
                           np.median(np.linalg.norm(H_A[:, wcols], axis=1))),
                       seconds=time.time() - t0)
            rows.append(row)
            out_path.write_text(json.dumps(rows, indent=2))
            print(f"[{held}/{arm}] nmae {row['nmae_pools']['wald89']:.5f}  "
                  f"prec89 {prec89:.5f}  fid_prec {row['fid_precision']:.5f}  "
                  f"jac89 {row['jac']['89']:.5f}  digest {row['digest']} "
                  f"({time.time() - t0:.0f}s)", flush=True)
            del final, H_A
        del op

    out_path.write_text(json.dumps(rows, indent=2))
    print(f"\n-> {out_path}")


if __name__ == "__main__":
    main()
