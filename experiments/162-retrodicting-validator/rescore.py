"""162 — second endpoint pass over the PLOLO folds.

`run_arms.py` recorded only SCALE-INVARIANT endpoints (`pds`, and `sse = 1 - cos_agg^2`
evaluated at the best global amplitude `gstar`).  Three of the six board members are not
scale-invariant, so that endpoint set structurally cannot price them:

  * `expr_mse_unbiased_capped_norm` is QUADRATIC in predicted amplitude (AGENTS.md 4.1).
    `sse` divides the amplitude out by construction, so an arm that only rescales the effect
    matrix reads `dsse = 0` and moves board `mse` a lot.  a19 is exactly that arm and it lost
    `-0.0959` scaled `mse`.  A validator built on `sse` is blind to the single largest `mse`
    move in the board record, by construction and not by accident.
  * `de_wilcoxon_lfc_nmae` is an absolute error on log2 fold changes.
  * `de_wilcoxon_direction_reach_raw` reads an ORDER, not a magnitude.

This pass adds, on the same cached folds and the same arm configs:

  `nmse_row`  mean_i ||P_i - T_i||^2 / ||T_i||^2      -- scale-SENSITIVE, per-row like the board
  `nmse_agg`  sum_i ||P_i - T_i||^2 / sum_i ||T_i||^2
  `nmae_k`    mean_i [ mean_{g in pool_i} |P-T| / mean_{g in pool_i} |T| ]
  `reach_k`   mean_i k*_i / K,  k*_i = max{k : purity of the first k >= 0.9}

`pool_i` is the top-K genes of row `i` by `|truth|` -- a stand-in for the reference's DE pool,
with K FIXED so that AGENTS.md 11.1 rule 9 (`n_real` in a denominator is unreadable on our
surfaces) is neutralised across arms rather than inherited.  K is a knob, so K in
{120, 240, 480} are all recorded and the retrodiction is required to be stable across them.

`reach_k` orders the pool by `|P|`.  The metric spec (BRIEF.md 5) orders by the submission's
adjusted p-value and reaches `|lfc|` only as the FOURTH tiebreak, so `reach_k` prices the
ordering's SIGN CONTENT and not the ordering the board actually walks.  That is a known,
stated substitution, and `retrodict.py` tests it rather than assuming it.

KNOWN-ANSWER GATE (AGENTS.md 11.1 rule 10): this script recomputes `pds`, `sse` and `gstar`
from scratch and requires that they reproduce `work/arms.json`.  That catches a fold cache,
operator, CD4 or config drift between the two passes, and it can fail: it compares two
independent process executions, not an algebraic identity.

The gate is a RELATIVE tolerance, not equality, and the reason is measured rather than
assumed.  Pass 1 and pass 2 disagree in the last few ulp on some cells (worst seen:
`k562_gwps/a2-prune` `cos` 0.031352451324800565 vs 0.03135245137817119, 1.7e-9 relative)
because `Operator.apply` runs a threaded BLAS `gemm` whose reduction order is not fixed across
processes.  Within one process it IS deterministic -- `run_arms.py`'s `a15` / `a15-dup` pair
is bit-identical on all 20 rows.  A configuration drift would show up at 1e-3 and larger, so a
1e-6 tolerance still fails loudly on everything the gate exists to catch, while 1e-9 float
noise does not stop the run.  The observed maximum is printed at the end; if it ever
approaches the tolerance that number, not the PASS, is the thing to read.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import plolo
from plolo import (FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver, endpoints,
                   solve_fraction)
from run_arms import CONFIGS, effects_for, frames_for

WORK = plolo.WORK
KS = (120, 240, 480)
PURITY = 0.9

# Gated at 1e-6.  Measured pass1-vs-pass2 worst case: `pds` 0.0 (a rank statistic, and it is
# bit-identical on all 165 rows), `sse` 4.6e-12, `cos` 2.6e-9, `gstar` 3.5e-9, `rr` 2.1e-7.
GATE_FIELDS = ("pds", "cos", "rr", "sse", "gstar")
# NOT gated, because they are genuinely irreproducible on ONE arm and the mechanism is known:
# `a23` runs `mass_consistent_effect`, a fixed-point iteration with an `rtol=1e-12` exit test.
# A 2e-9 shift in its input changes which iteration exits, which flips the sign of genes whose
# delivered effect is ~0.  Measured drift: `sgn` 6.1e-3, `sgn_w` 8.7e-3 on a23 and EXACTLY 0.0
# on every other arm.  So any sign-keyed statistic read off a mass-fixed prediction -- which
# includes `reach_k` below, and board `fid`/`reach`/`jac` -- carries a ~1% reproducibility
# floor on that arm that has nothing to do with the model.
SIGN_FIELDS = ("sgn", "sgn_w")


def de_endpoints(P, T, ks=KS, purity=PURITY):
    """`nmae_k` and `reach_k` over a fixed top-K-by-|truth| pool, per row then averaged."""
    out = {}
    absT = np.abs(T)
    order = np.argsort(-absT, axis=1)          # genes by |truth| descending, per row
    sgn_ok_all = np.sign(P) == np.sign(T)
    for K in ks:
        k = min(K, T.shape[1])
        idx = order[:, :k]
        rows = np.arange(T.shape[0])[:, None]
        t = T[rows, idx]
        p = P[rows, idx]
        denom = np.abs(t).mean(1)
        nmae = np.mean(np.abs(p - t).mean(1) / np.maximum(denom, 1e-300))

        # reach: re-order the pool by |P| descending, walk, deepest prefix with purity >= 0.9
        sub = np.argsort(-np.abs(p), axis=1)
        ok = sgn_ok_all[rows, idx][rows, sub].astype(np.float64)
        run = np.cumsum(ok, axis=1) / np.arange(1, k + 1)[None, :]
        clears = run >= purity
        depth = np.where(clears.any(1), k - np.argmax(clears[:, ::-1], axis=1), 0)
        out[f"nmae_k{K}"] = float(nmae)
        out[f"reach_k{K}"] = float(np.mean(depth / k))
    return out


def mse_endpoints(P, T):
    """Scale-SENSITIVE normalised squared error, per-row and aggregate."""
    num = ((P - T) ** 2).sum(1)
    den = (T * T).sum(1)
    return dict(nmse_row=float(np.mean(num / np.maximum(den, 1e-300))),
                nmse_agg=float(num.sum() / max(den.sum(), 1e-300)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--arms", nargs="*", default=[a for a in CONFIGS if a != "a15-dup"])
    ap.add_argument("--out", default=str(WORK / "arms2.json"))
    ap.add_argument("--gate-rtol", type=float, default=1e-6)
    args = ap.parse_args()

    ref = {(r["fold"], r["arm"], r["space"]): r
           for r in json.loads((WORK / "arms.json").read_text())}
    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"], r["space"]) for r in rows}
    gate_checked = 0
    gate_worst = (0.0, "")
    sign_worst = (0.0, "")

    for held in args.folds:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        cd4 = cd4_for(fold["targets"])
        frames = frames_for(fold)
        ops: dict[str, Operator] = {}
        for arm in args.arms:
            cfg = CONFIGS[arm]
            if (held, arm, "chalA/delivered") in done:
                continue
            if cfg["op"] != "union_pruned" and FOLDS[held]["op_kind"] != "production":
                continue
            key = cfg["op"] if FOLDS[held]["op_kind"] == "production" else "held"
            if key not in ops:
                path = (PROD_OPS / f"whitener_{cfg['op']}.npz"
                        if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
                ops[key] = Operator(path)
            final, caud = effects_for(fold, cfg, ops[key], cd4)

            for space, mat in [("intended", final)] + [
                    (f"{n}/delivered", deliver(final, b, g, cfg["mass"]))
                    for n, b, g in frames]:
                P = mat[:, fold["scored"]]
                e = endpoints(P, truth)
                # --- known-answer gate against pass 1 ---
                prev = ref.get((held, arm, space))
                if prev is not None:
                    for k in GATE_FIELDS:
                        rel = abs(prev[k] - e[k]) / max(abs(prev[k]), 1e-12)
                        if rel > gate_worst[0]:
                            gate_worst = (rel, f"{held}/{arm}/{space}/{k}")
                        if rel > args.gate_rtol:
                            raise SystemExit(
                                f"KA GATE FAIL {held}/{arm}/{space}/{k}: "
                                f"pass1 {prev[k]!r} != pass2 {e[k]!r} (rel {rel:.2e})")
                    for k in SIGN_FIELDS:
                        rel = abs(prev[k] - e[k]) / max(abs(prev[k]), 1e-12)
                        if rel > sign_worst[0]:
                            sign_worst = (rel, f"{held}/{arm}/{space}/{k}")
                    gate_checked += 1
                rows.append(dict(fold=held, arm=arm, space=space, **e,
                                 **mse_endpoints(P, truth), **de_endpoints(P, truth),
                                 **caud, n_targets=len(fold["targets"]),
                                 n_scored=int(len(fold["scored"]))))
                r = rows[-1]
                print(f"  {held:<14} {arm:<9} {space:<18} pds {r['pds']:.6f} "
                      f"nmse {r['nmse_row']:.5f} nmae240 {r['nmae_k240']:.5f} "
                      f"reach240 {r['reach_k240']:.5f}", flush=True)
            out_path.write_text(json.dumps(rows, indent=1))
        del ops
    print(f"KA gate: {gate_checked} endpoint rows reproduced pass 1 within rtol "
          f"{args.gate_rtol:.0e}; worst {gate_worst[0]:.2e} at {gate_worst[1]}")
    print(f"sign-field drift (not gated): worst {sign_worst[0]:.2e} at {sign_worst[1]}")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
