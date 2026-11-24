"""216 Phase 2.1 — take the CD4 dose to its limit.

207 stopped reporting at `f = 0.80` (its cosine argmax) and shipped 0.34.  The mandate asks
for the curve well past that, because the board's two leaders moved their cosine by x2.06 in
five days and the CD4 dose is the only lever we hold that has ever moved cosine.

Differences from `sweep207.py`, all deliberate:

  * a FINER grid through the turnover (0.70 -> 2.0 in eleven steps) and a much longer tail
    (to f = 40), so "where does it turn over" is answered by measurement and not by the two
    points either side of it;
  * both cosine functionals are carried at both tags.  207 published `cos_agg` on the a23
    (mass-fixed) tag; 204 and 212 published the macro `cos` on the a15 tag.  Phase 0 showed
    that comparing across tags is what makes a known-answer row look like it failed;
  * ONE operator resident at a time (207 held four; that is where its 4.6 GB went and where
    this experiment's ceiling would have gone);
  * the target-shuffled falsifier is run at EVERY dose in the extended range, not only at
    the twelve 207 reported.

    uv run python sweep216.py
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
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (E062, FOLDS, OFFICIAL, PROD_OPS, Operator, build_fold,  # noqa: E402
                   cd4_for, deliver, endpoints, solve_fraction)

ALPHA = 0.982
GAMMA = 0.65
LAM = 0.0

# The grid.  207's points are all present so every one of its rows is a known answer here.
DOSE_F = [0.0, 0.05, 0.10, 0.14, 0.20, 0.28, 0.34, 0.40, 0.47, 0.55, 0.62, 0.70, 0.80,
          0.90, 1.00, 1.10, 1.20, 1.30, 1.45, 1.60, 1.80, 2.00, 2.50, 3.00, 4.00, 5.50,
          7.00, 10.0, 15.0, 25.0, 40.0]
SHUF_F = [0.05, 0.10, 0.20, 0.34, 0.40, 0.55, 0.80, 1.00, 1.30, 2.00, 3.00, 7.00, 15.0]


def shuffled_cd4(cd4, seed: int = 20726):
    """207's falsifier-1 construction, imported by value so the seed and the derangement
    rule are identical and its twelve published rows are reproducible here."""
    rng = np.random.default_rng(seed)
    n = cd4.mean_effect.shape[0]
    carry = np.flatnonzero(np.asarray(cd4.target_condition_count) == 3)
    perm = np.arange(n)
    shuffled = carry.copy()
    while True:
        rng.shuffle(shuffled)
        if len(carry) < 2 or not np.any(shuffled == carry):
            break
    perm[carry] = shuffled
    return E062.CD4Prepared(
        mode=cd4.mode + "+shuffled",
        mean_effect=np.asarray(cd4.mean_effect)[perm].copy(),
        gene_votes=np.asarray(cd4.gene_votes)[perm].copy(),
        target_condition_count=np.asarray(cd4.target_condition_count)[perm].copy(),
        condition_details=cd4.condition_details)


def _per_target(P, T):
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    pds = ((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1)
    return pds, own


def score_effects(final, fold, truth):
    out, per, perc = {}, {"a15": [], "a23": []}, {"a15": [], "a23": []}
    for tag, mass in (("a23", "minnorm"), ("a15", None)):
        acc = {}
        for ctx in ("A", "B", "C"):
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            gate = (cpm > plolo.MIN_CPM) & (b > 0)
            dl = deliver(final, b, gate, mass)
            e = endpoints(dl[:, fold["scored"]], truth)
            e["cos_agg"] = float(np.sqrt(max(0.0, 1.0 - e["sse"])))
            acc[ctx] = e
            p, c = _per_target(dl[:, fold["scored"]], truth)
            per[tag].append(p)
            perc[tag].append(c)
        for k in acc["A"]:
            out[f"{tag}_{k}"] = float(np.mean([acc[c][k] for c in "ABC"]))
    out.update({k: out[f"a23_{k}"] for k in ("pds", "cos", "cos_agg", "sgn", "sgn_w", "sse")})
    out["_pds_pt"] = np.mean(per["a23"], axis=0).tolist()
    out["_cos_pt"] = np.mean(perc["a23"], axis=0).tolist()
    out["_cos_pt_a15"] = np.mean(perc["a15"], axis=0).tolist()
    out["median_row_norm_final"] = float(np.median(np.linalg.norm(final, axis=1)))
    return out


def main() -> None:
    t0 = time.time()
    out_path = WORK / "sweep216.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    arms = [("f%g" % v, "cd4", v) for v in DOSE_F]
    arms += [("shuf%g" % v, "shuffled", v) for v in SHUF_F]

    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        shuf = shuffled_cd4(cd4)
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        whitened = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)
        del op
        gc.collect()

        for name, kind, f_val in arms:
            if (held, name) in done:
                continue
            t1 = time.time()
            resid = cd4 if kind == "cd4" else shuf
            if f_val > 0:
                frac, per_unit = solve_fraction(whitened, resid, f_val)
                final, preserve = E062.add_post_whiten_residual(
                    whitened, resid, frac, fold["wcols"])
                ps = float(preserve["post_add_global_preserve_scale"])
            else:
                frac, per_unit, ps = 0.0, None, 1.0
                final = whitened.copy()
            final[:, fold["panel_positions"]] = 0.0
            e = score_effects(final, fold, truth)
            rows.append(dict(fold=held, arm=name, kind=kind, f=f_val,
                             n_targets=len(fold["targets"]), fraction=float(frac),
                             f_per_unit=(float(per_unit) if per_unit else None),
                             preserve_scale=ps, seconds=round(time.time() - t1, 1), **e))
            out_path.write_text(json.dumps(rows))
            print(f"[{held}/{name}] frac {frac:.5f} pds {e['pds']:.6f} "
                  f"cos_agg {e['cos_agg']:.6f} cos15 {e['a15_cos']:.7f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
            del final
        del whitened, pre, cd4, shuf, fold, truth
        gc.collect()
    print(f"done in {time.time()-t0:.0f}s, {len(rows)} rows")


if __name__ == "__main__":
    main()
