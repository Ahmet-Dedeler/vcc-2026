"""236 — THE CONTROL 233 DID NOT RUN, and the one this stack most needs.

The four October sources carry almost the SAME rows (measured: 28/28/27/25 singles, union 30,
on `pess/orion_hct116`).  So the whole stack gain could be nothing but **more residual mass on
a favourable 30-row subset** -- a row-selection effect, not new biology.  218 built exactly the
control that separates those two: `compose218.restrict_rows(cd4, stack)` -- **CD4's own effect,
on the stack's own rows, at the stack's own dose, sized the same carrying-row way.**  If the
gain is "more residual on these rows", this reproduces it and the October sources are worth
nothing beyond a dose change we already control.

Also run: `rand` -- CD4's effect on a RANDOM row set of the same size, which separates "these
particular rows" from "any rows at all".

    uv run python rowctl236.py
"""
from __future__ import annotations

import argparse
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
R233 = ROOT / "experiments/233-october-retention"
sys.path.insert(0, str(R233))

import panels233 as P  # noqa: E402
import score233 as S  # noqa: E402
import fold233  # noqa: E402
import plolo  # noqa: E402
from plolo import E062, Operator  # noqa: E402

from stack236 import FORBIDDEN, OCTOBER_ONLY, DOSES  # noqa: E402

RAND_SEED = 2361


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="*", default=["pess", "real"])
    ap.add_argument("--folds", nargs="*", default=list(P.CANCER_FOLDS))
    ap.add_argument("--out", default=str(WORK / "rowctl236.json"))
    args = ap.parse_args()

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["panel"], r["fold"], r["arm"]) for r in rows}
    resid218, compose218 = S._resid_module()

    ctrl_mod = __import__("ctrl228")
    ctrl = {c: ctrl_mod.control_block(c) for c in S.CONTEXTS}
    bctx = {c: plolo.challenge_bctrl(c) for c in S.CONTEXTS}
    gate = {}
    for c in S.CONTEXTS:
        e = np.expm1(bctx[c])
        gate[c] = (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (bctx[c] > 0)

    for panel in args.panels:
        S.bind_panel(panel, resid218)
        for held in args.folds:
            t0 = time.time()
            fold = fold233.build(panel, held)
            # FALSIFIER 5: the operator is a function of the panel and fails silently
            print(f"[{panel}/{held}] operator {fold233.verify_op(panel, held)}", flush=True)
            sc = np.asarray(fold["scored"])
            wcols = np.asarray(fold["wcols"])
            Ts = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
            targets = fold["targets"]
            cd4 = S.cd4_prepared(panel, targets)
            h1 = compose218.union_rows(
                [resid218.source_residual(n, targets, "gene_norm") for n in S.UNION])
            op = Operator(fold["op_path"])
            pre = fold233.consensus(fold, pm=True, only=None)[0]
            whitened = op.apply(pre * S.ALPHA, wcols, S.GAMMA, S.LAM)
            frac_cd4, _ = plolo.solve_fraction(whitened, cd4, S.CD4_F)
            frac_u, _ = compose218.solve_fraction_carrying(whitened, h1, S.UNION_G)
            base_resid = compose218.merged(cd4, frac_cd4, h1, frac_u)

            from mapfit import to_lfc
            from phase23 import POOL_N, topn
            cpm_sc = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(
                np.expm1(fold["b_ctrl"]).sum())
            LR = to_lfc(np.asarray(fold["truth"], dtype=np.float64),
                        fold["b_ctrl"], gate["A"])[:, sc]
            pools = {"wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm_sc, 0)), POOL_N),
                     "lfc89": topn(np.abs(LR), POOL_N), "h89": topn(np.abs(Ts), POOL_N)}
            real_de = {n: topn(np.abs(Ts), n) for n in S.N_REALS}

            names = tuple(n for n in OCTOBER_ONLY if n not in FORBIDDEN.get(held, set()))
            stack = compose218.union_rows(
                [resid218.source_residual(n, targets, "gene_norm") for n in names])
            rowcd4 = compose218.restrict_rows(cd4, stack)
            # a random row set of the same size, carrying CD4's own effect
            live = np.asarray(stack.target_condition_count) > 0
            rng = np.random.default_rng(RAND_SEED)
            pick = rng.choice(len(live), size=int(live.sum()), replace=False)
            fake = E062.CD4Prepared(
                mode="rand_rows", mean_effect=np.zeros_like(np.asarray(stack.mean_effect)),
                gene_votes=np.zeros(np.asarray(stack.gene_votes).shape, dtype=np.uint8),
                target_condition_count=np.zeros(live.shape, dtype=np.uint8),
                condition_details=[])
            np.asarray(fake.target_condition_count)[pick] = 1
            randcd4 = compose218.restrict_rows(cd4, fake)

            arms = {"stack": stack, "rowcd4": rowcd4, "randcd4": randcd4}
            for arm_name, res in arms.items():
                carrying = int((np.asarray(res.target_condition_count) > 0).sum())
                for g in DOSES:
                    key = f"{arm_name}@{g:g}"
                    if (panel, held, key) in done:
                        continue
                    fx, per_unit = compose218.solve_fraction_carrying(whitened, res, g)
                    add = compose218.merged(base_resid, 1.0, res, fx)
                    final, _ = E062.add_post_whiten_residual(whitened, add, 1.0, wcols)
                    final[:, fold["panel_positions"]] = 0.0
                    m = S.score_arm(final, fold, bctx, gate, ctrl, Ts, LR, pools, real_de,
                                    np.random.default_rng(7), None, full=False)
                    rows.append(dict(panel=panel, fold=held, arm=key, dose=g,
                                     members=list(names), carrying_rows=carrying,
                                     per_unit_carrying=per_unit, fraction=float(fx),
                                     seconds=round(time.time() - t0, 1), **m))
                    out_path.write_text(json.dumps(rows, indent=2))
                    print(f"[{panel}/{held}/{key}] cos {m['cos_agg']:.6f} "
                          f"pds {m['pds']:.6f} carry {carrying}", flush=True)
                    del final
            del op, whitened, pre, stack, rowcd4, randcd4
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
