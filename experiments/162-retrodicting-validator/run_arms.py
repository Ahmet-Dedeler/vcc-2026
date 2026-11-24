"""162 — run every reconstructible past arm through PLOLO and record its endpoints.

One row per (fold, arm-config).  Deltas are formed afterwards in `retrodict.py`, so this
script never sees a board number and cannot be tuned against one.

KNOWN-ANSWER ROWS (AGENTS.md 11.1 rule 10 — a screen with no row whose answer is already
known cannot detect a whole class of silent bug):

  1. `a15` re-run under a second config object that is field-identical must reproduce every
     endpoint bit-for-bit.  Catches accidental state carried between arms.
  2. `identity-intended` vs `identity-delivered`: 156 measured +0.00324/+0.00382 raw pds on
     its own two surfaces.  Sign and rough size must reproduce here or the renderer is
     wired wrong.
  3. The four emission arms (a8, a10, a21, a22) are not run at all: they cannot change the
     effect matrix, and the board confirms `pds`/`mse` moved by exactly 0.0 on a21/a22 to
     sixteen significant digits.  They enter the retrodiction table as exact zeros.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import plolo
from plolo import (FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver, endpoints,
                   solve_fraction)

WORK = plolo.WORK

# alpha / (gamma, lam) / cd4 / mass_fix / operator, exactly as the shipped meta.json files record
CONFIGS: dict[str, dict] = {
    # v11-era: full union operator, pre-096 amplitude
    "v11":        dict(alpha=1.0796, gamma=0.65, lam=0.0, cd4_frac=0.05,  mass=None, op="union"),
    "a2-prune":   dict(alpha=1.0796, gamma=0.65, lam=0.0, cd4_frac=0.05,  mass=None, op="union_pruned"),
    # a10 == a8 == a2-prune in geometry (emission-only changes)
    "a10":        dict(alpha=1.0796, gamma=0.65, lam=0.0, cd4_frac=0.05,  mass=None, op="union_pruned"),
    "a11":        dict(alpha=0.838,  gamma=0.65, lam=0.0, cd4_frac=0.05,  mass=None, op="union_pruned"),
    "a12":        dict(alpha=1.0108, gamma=0.65, lam=0.0, cd4_frac=0.05,  mass=None, op="union_pruned"),
    "a14":        dict(alpha=0.982,  gamma=0.65, lam=0.0, cd4_frac=0.05,  mass=None, op="union_pruned"),
    "a15":        dict(alpha=0.982,  gamma=0.65, lam=0.0, cd4_f=0.10,     mass=None, op="union_pruned"),
    "a19":        dict(alpha=0.982,  gamma=0.80, lam=0.007, cd4_f=0.10,   mass=None, op="union_pruned"),
    "a23":        dict(alpha=0.982,  gamma=0.65, lam=0.0, cd4_f=0.10,     mass="minnorm", op="union_pruned"),
    # known-answer duplicate of a15
    "a15-dup":    dict(alpha=0.982,  gamma=0.65, lam=0.0, cd4_f=0.10,     mass=None, op="union_pruned"),
}


def effects_for(fold, cfg, op, cd4):
    pre = np.asarray(fold["pre"], dtype=np.float64) * cfg["alpha"]
    whitened = op.apply(pre, fold["wcols"], cfg["gamma"], cfg["lam"])
    if "cd4_f" in cfg:
        frac, per_unit = solve_fraction(whitened, cd4, cfg["cd4_f"])
    else:
        frac, per_unit = cfg["cd4_frac"], None
    final, _ = plolo.E062.add_post_whiten_residual(whitened, cd4, frac, fold["wcols"])
    final[:, fold["panel_positions"]] = 0.0
    eff_f = frac * (per_unit if per_unit else
                    float(np.median(np.linalg.norm(
                        np.where(cd4.gene_votes > 0, cd4.mean_effect, 0.0), axis=1)))
                    / (float(np.median(np.linalg.norm(whitened, axis=1))) or 1.0))
    return final, dict(cd4_fraction=float(frac), effective_f=float(eff_f))


def frames_for(fold):
    """The two render frames the validator is allowed to use, and why both are tried.

    `library`   : the held library's own control profile.  Self-consistent -- prediction
                  and truth live in the same frame -- but its render residual is NOT the
                  residual the board sees, because the frame is a different cell line on
                  a different gene support.
    `challenge` : contexts A/B/C, i.e. the exact frames production renders into.  The
                  render residual is then the real one; the truth is still the held
                  library's, which is fine because the effect matrix is context-independent
                  by construction and the truth is a log-fold-change.
    """
    out = [("library", fold["b_ctrl"], fold["gate"])]
    for ctx in ("A", "B", "C"):
        b = plolo.challenge_bctrl(ctx)
        expr = np.expm1(b)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        out.append((f"chal{ctx}", b, (cpm > plolo.MIN_CPM) & (b > 0)))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--arms", nargs="*", default=list(CONFIGS))
    ap.add_argument("--out", default=str(WORK / "arms.json"))
    args = ap.parse_args()

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"], r["space"]) for r in rows}

    for held in args.folds:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        cd4 = cd4_for(fold["targets"])
        frames = frames_for(fold)
        ops: dict[str, Operator] = {}
        for arm in args.arms:
            cfg = CONFIGS[arm]
            if (held, arm, "library/delivered") in done:
                continue
            if cfg["op"] != "union_pruned" and FOLDS[held]["op_kind"] != "production":
                continue  # only the leak-free production fold can express an operator swap
            key = cfg["op"] if FOLDS[held]["op_kind"] == "production" else "held"
            if key not in ops:
                path = (PROD_OPS / f"whitener_{cfg['op']}.npz"
                        if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
                ops[key] = Operator(path)
                print(f"[{held}] loaded operator {path.name} "
                      f"maxeig {ops[key].evals.max():.4f}", flush=True)
            final, caud = effects_for(fold, cfg, ops[key], cd4)
            e = endpoints(final[:, fold["scored"]], truth)
            rows.append(dict(fold=held, arm=arm, space="intended", **e, **caud,
                             n_targets=len(fold["targets"]),
                             n_scored=int(len(fold["scored"]))))
            print(f"  {held:<14} {arm:<9} {'intended':<18} pds {e['pds']:.6f} "
                  f"cos {e['cos']:.5f} sse {e['sse']:.6f} sgnw {e['sgn_w']:.5f}", flush=True)
            for fname, b, gate in frames:
                diag = plolo.render_diagnostics(final, b, gate)
                mat = deliver(final, b, gate, cfg["mass"])
                e = endpoints(mat[:, fold["scored"]], truth)
                space = f"{fname}/delivered"
                rows.append(dict(fold=held, arm=arm, space=space, **e, **caud, **diag,
                                 n_targets=len(fold["targets"]),
                                 n_scored=int(len(fold["scored"]))))
                print(f"  {held:<14} {arm:<9} {space:<18} pds {e['pds']:.6f} "
                      f"cos {e['cos']:.5f} sse {e['sse']:.6f} sgnw {e['sgn_w']:.5f} "
                      f"| massratio {diag['mass_ratio']:.5f} "
                      f"resid {diag['render_residual_pct']:.2f}%", flush=True)
            out_path.write_text(json.dumps(rows, indent=1))
        del ops
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
