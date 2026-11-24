"""207 Phase 2 — the FORM: CD4 as a FOURTH DONOR LIBRARY vs the post-whitening residual.

GISL ships the Marson/Zhu CD4 atlas as a **library** (a fourth member of the consensus,
whitened with the rest) and published a board-measured +0.01325 overall for it.  We ship it
as a **post-whitening residual**.  Those are different estimators, not two doses of one.

Construction of the library arm, faithful to the incumbent's own conventions so that the arm
is ONE variable (rule 10):

    build_baseline sums the sources that carry a target and divides by the NUMBER of them,
    counting a zero on a gene the source does not measure (039's implicit coverage shrinkage,
    deliberately preserved).  So the fourth-library form is

        pre4 = (pre * votes + w * cd4_effect) / (votes + w)      on rows CD4 carries
        pre4 = pre                                               elsewhere

    `votes` is recomputed from the donors' own `perts.txt`, and `pre * votes` reconstructs the
    numerator exactly because `pre` is that numerator divided by `votes`.  KNOWN-ANSWER ROW:
    `w = 0` must reproduce `f0`'s endpoints to machine precision.

`cd4.mean_effect` is already on the donors' EFFECT_SCALE convention (062 `prepare_cd4_mean`
rescales each condition to EFFECT_SCALE exactly as `prepare_main_source` does), so `w = 1` is
a genuine equal-weight fourth member.

    uv run python library207.py
"""
from __future__ import annotations

import argparse
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
                   solve_fraction)
from sweep207 import ALPHA, GAMMA, LAM, sha8, score_effects  # noqa: E402

LIB_W = [0.0, 0.0625, 0.125, 0.1875, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875,
         1.0, 1.25, 1.5, 2.0, 3.0]
COMPOSE = [(w, f) for w in (0.25, 0.375, 0.5, 0.625, 0.75, 1.0)
           for f in (0.05, 0.10, 0.15, 0.20, 0.28, 0.40, 0.55)]


def votes_for(held: str, targets) -> np.ndarray:
    v = np.zeros(len(targets), dtype=np.float64)
    for name in FOLDS[held]["sources"]:
        perts = set((plolo.DONOR_PATHS[name] / "perts.txt").read_text().split())
        v += np.array([t in perts for t in targets], dtype=np.float64)
    return v


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--out", default=str(WORK / "library207.json"))
    args = ap.parse_args()

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in args.folds:
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
        op_meta = dict(path=str(op_path.relative_to(ROOT)), sha8=sha8(op_path))
        n0 = float(np.median(np.linalg.norm(pre[:, fold["wcols"]], axis=1)))

        arms = [(f"lib{w:g}", w, 0.0) for w in LIB_W]
        arms += [(f"lib{w:g}+res{f:g}", w, f) for w, f in COMPOSE]
        for arm, w, resf in arms:
            if (held, arm) in done:
                continue
            t0 = time.time()
            rowsel = carry & (votes > 0)
            pre4 = pre.copy()
            if w > 0:
                pre4[rowsel] = ((pre[rowsel] * votes[rowsel, None] + w * cd4e[rowsel])
                                / (votes[rowsel, None] + w))
            n1 = float(np.median(np.linalg.norm(pre4[:, fold["wcols"]], axis=1)))
            # `alpha*`: match the incumbent's PRE-whitening median row norm, so the arm is a
            # direction change and not an amplitude change (203's re-solve, applied up front).
            alpha_star = ALPHA * n0 / n1 if n1 else ALPHA
            for tag, a in (("", ALPHA), ("@a*", alpha_star)):
                name = arm + tag
                if (held, name) in done:
                    continue
                wh = op.apply(pre4 * a, fold["wcols"], GAMMA, LAM)
                if resf > 0:
                    frac, _ = solve_fraction(wh, cd4, resf)
                    final, _ = E062.add_post_whiten_residual(wh, cd4, frac, fold["wcols"])
                else:
                    frac, final = 0.0, wh
                final[:, fold["panel_positions"]] = 0.0
                e = score_effects(final, fold, truth)
                aud = dict(lib_w=w, residual_f=resf, alpha=a, alpha_star=alpha_star,
                           pre_norm_ratio=n1 / n0, cd4_rows=int(rowsel.sum()),
                           rows_votes0=int((votes == 0).sum()), fraction=float(frac),
                           operator=op_meta)
                rows.append(dict(fold=held, arm=name, n_targets=len(fold["targets"]),
                                 note=f"CD4 as fourth library member w={w:g}, "
                                      f"post-whiten residual f={resf:g}, alpha={a:.4f}",
                                 **e, audit=aud, seconds=round(time.time() - t0, 1)))
                out_path.write_text(json.dumps(rows, indent=2))
                print(f"[{held}/{name}] a {a:.4f} pds {e['pds']:.6f} "
                      f"cos_agg {e['cos_agg']:.6f} cos {e['cos']:.6f} "
                      f"q_w {e['a15_sgn_w']:.6f} prenorm {n1/n0:.4f}", flush=True)
        del op


if __name__ == "__main__":
    main()
