"""216 Phase 2.2/2.3 — re-open the LIBRARY form against 207's Phase 2.

207 measured the library form losing on all 57 settings, with `nmae` exploding because the
**delivered** median row norm collapsed to **0.523x** a15's even at its re-solved alpha.  The
re-solve it ran matched the **pre-whitening** norm (`library207.py:103`), and 207 §5.4 point 2
says in as many words that this is not the norm that collapsed: "the delivered median row norm
is still 0.523x a15's, because averaging in a fourth member shrinks rows non-uniformly and the
renderer is nonlinear".

So the honest re-test is an alpha solved against the **DELIVERED** norm -- the quantity the
board actually reads -- and that is the third alpha policy here:

    `@a0`    alpha = 0.982, 207's baseline row (known-answer: must reproduce `library207.json`)
    `@apre`  207's pre-whitening norm match     (known-answer: same)
    `@adel`  NEW: alpha solved so the delivered median row norm equals the incumbent's

`lib1` is the mandate's 1/k route: `build_baseline` divides by the number of carrying sources,
so `w = 1` makes CD4 one more equal member -- a genuine FOURTH equal line on `kolf2` (three
sources) and a third on the other folds.  159's law predicts only sqrt(4/3) = 1.155x from
3 -> 4 lines, which is reported alongside the measurement.

    uv run python form216.py
"""
from __future__ import annotations

import gc
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   deliver, solve_fraction)
from sweep216 import ALPHA, GAMMA, LAM, score_effects  # noqa: E402

LIB_W = [0.0, 0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
COMPOSE = [(w, f) for w in (0.25, 0.5, 1.0) for f in (0.10, 0.20, 0.34, 0.55)]


def votes_for(held: str, targets) -> np.ndarray:
    v = np.zeros(len(targets), dtype=np.float64)
    for name in FOLDS[held]["sources"]:
        perts = set((plolo.DONOR_PATHS[name] / "perts.txt").read_text().split())
        v += np.array([t in perts for t in targets], dtype=np.float64)
    return v


def delivered_norm(final, fold) -> float:
    """Median row norm of the object the SCORER reads, averaged over the three frames.

    This is the quantity 207 §5.4 measured at 0.523x for the library form; matching it is
    what the `@adel` policy does.  `minnorm` is production's mass fix (the a23 tag).
    """
    vals = []
    for ctx in ("A", "B", "C"):
        b = plolo.challenge_bctrl(ctx)
        expr = np.expm1(b)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        gate = (cpm > plolo.MIN_CPM) & (b > 0)
        dl = deliver(final, b, gate, "minnorm")
        vals.append(float(np.median(np.linalg.norm(dl[:, fold["scored"]], axis=1))))
    return float(np.mean(vals))


def build_final(pre4, alpha, op, fold, cd4, resf):
    wh = op.apply(pre4 * alpha, fold["wcols"], GAMMA, LAM)
    if resf > 0:
        frac, _ = solve_fraction(wh, cd4, resf)
        final, _ = E062.add_post_whiten_residual(wh, cd4, frac, fold["wcols"])
    else:
        frac, final = 0.0, wh
    final[:, fold["panel_positions"]] = 0.0
    return final, float(frac)


def solve_alpha_delivered(pre4, op, fold, cd4, resf, target, a0=ALPHA, tol=2e-4, it=12):
    """Secant solve on log(delivered median row norm) -- monotone increasing in alpha."""
    xs = [a0, a0 * 1.15]
    ys = []
    for a in xs:
        f, _ = build_final(pre4, a, op, fold, cd4, resf)
        ys.append(np.log(delivered_norm(f, fold)) - np.log(target))
        del f
    for _ in range(it):
        if abs(ys[-1]) < tol:
            break
        denom = (ys[-1] - ys[-2])
        if denom == 0:
            break
        a = xs[-1] - ys[-1] * (xs[-1] - xs[-2]) / denom
        a = float(np.clip(a, 0.2 * a0, 6.0 * a0))
        f, _ = build_final(pre4, a, op, fold, cd4, resf)
        xs.append(a)
        ys.append(np.log(delivered_norm(f, fold)) - np.log(target))
        del f
    return float(xs[-1]), float(np.exp(ys[-1]))


def main() -> None:
    t0 = time.time()
    out_path = WORK / "form216.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        votes = votes_for(held, fold["targets"])
        carry = np.asarray(cd4.target_condition_count) == 3
        cd4e = np.asarray(cd4.mean_effect, dtype=np.float64)
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        n0 = float(np.median(np.linalg.norm(pre[:, fold["wcols"]], axis=1)))

        # the incumbent: the shipped residual at f = 0.10, alpha 0.982.  Its delivered
        # median row norm is the target every `@adel` arm is solved onto.
        inc_final, _ = build_final(pre, ALPHA, op, fold, cd4, 0.10)
        target_norm = delivered_norm(inc_final, fold)
        inc = score_effects(inc_final, fold, truth)
        del inc_final
        gc.collect()
        print(f"[{held}] incumbent f=0.10: delivered median row norm {target_norm:.6f} "
              f"cos_agg {inc['cos_agg']:.6f} pds {inc['pds']:.6f}", flush=True)

        arms = [(f"lib{w:g}", w, 0.0) for w in LIB_W]
        arms += [(f"lib{w:g}+res{f:g}", w, f) for w, f in COMPOSE]

        for arm, w, resf in arms:
            rowsel = carry & (votes > 0)
            pre4 = pre.copy()
            if w > 0:
                pre4[rowsel] = ((pre[rowsel] * votes[rowsel, None] + w * cd4e[rowsel])
                                / (votes[rowsel, None] + w))
            n1 = float(np.median(np.linalg.norm(pre4[:, fold["wcols"]], axis=1)))
            a_pre = ALPHA * n0 / n1 if n1 else ALPHA
            policies = [("@a0", ALPHA, None), ("@apre", a_pre, None), ("@adel", None, "solve")]
            for tag, a, mode in policies:
                name = arm + tag
                if (held, name) in done:
                    continue
                t1 = time.time()
                solve_resid = None
                if mode == "solve":
                    a, solve_resid = solve_alpha_delivered(
                        pre4, op, fold, cd4, resf, target_norm, a0=a_pre)
                final, frac = build_final(pre4, a, op, fold, cd4, resf)
                dn = delivered_norm(final, fold)
                e = score_effects(final, fold, truth)
                rows.append(dict(
                    fold=held, arm=name, lib_w=w, residual_f=resf, alpha=float(a),
                    alpha_pre=float(a_pre), alpha_policy=tag, fraction=frac,
                    pre_norm_ratio=float(n1 / n0),
                    delivered_norm=dn, delivered_ratio=float(dn / target_norm),
                    delivered_target=target_norm, solve_residual=solve_resid,
                    cd4_rows=int(rowsel.sum()), n_sources=len(FOLDS[held]["sources"]),
                    seconds=round(time.time() - t1, 1),
                    **{k: v for k, v in e.items() if not k.startswith("_")},
                    _cos_pt=e["_cos_pt"], _pds_pt=e["_pds_pt"]))
                out_path.write_text(json.dumps(rows))
                print(f"[{held}/{name}] a {a:.4f} delivered {dn/target_norm:.4f}x "
                      f"pds {e['pds']:.6f} cos_agg {e['cos_agg']:.6f} "
                      f"({time.time()-t0:.0f}s)", flush=True)
                del final
            del pre4
            gc.collect()
        del op, pre, cd4, fold, truth
        gc.collect()
    print(f"done in {time.time()-t0:.0f}s, {len(rows)} rows")


if __name__ == "__main__":
    main()
