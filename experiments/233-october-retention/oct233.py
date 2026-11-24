"""233 Phase 3 — price the FOUR SOURCES THAT ONLY EXIST IN OCTOBER.

`nadig_hepg2`, `nadig_jurkat`, `rpe1`, `k562_essential` are common-essential screens.  They
carry **0 of today's 300** targets, so PLOLO cannot price them at all today -- 225 called their
October value *"the single largest untested claim in this file."*  On a disjoint panel they
carry 30-48 each (measured, `census233.py`), and become priceable for the first time.

FORM.  218's, unchanged: an un-whitened post-whitening residual added on top of the shipped
stack (CD4 @ 0.34 + H1 union @ 0.28), at a dose `g`, scored on PLOLO after the renderer and
the scorer's read-back.  Primary endpoint the energy-weighted aggregate cosine; `pds` second.

FALSIFIER 4 — THE SIZING RULE.  These sources cover 10-16% of a fold's targets.  Production's
`solve_cd4_fraction` takes the median row norm over EVERY row, so below 50% coverage it
returns `f_per_unit = 0.0` and the residual is not small, it is IDENTICALLY ABSENT (223's
silent zero; the harness mirror at `162.../plolo.py:363-369`).  A disjoint panel is exactly
what produces low-coverage sources, so every dose here is solved on the CARRYING rows
(`compose218.solve_fraction_carrying`) and the all-rows solve is reported alongside to show
it is zero.

LEAK AUDIT.  `k562_essential` is a Replogle K562 screen and `k562_gwps` is the same cell line,
so it is NEVER used on the fold that holds K562 out as truth.  Every other pairing is a
different line.  A target-deranged copy of each source is run at the same doses as the null.

    uv run python oct233.py --panels real pess
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
sys.path.insert(0, str(HERE))

import panels233 as P  # noqa: E402
import score233 as S  # noqa: E402
import fold233  # noqa: E402
import plolo  # noqa: E402
from plolo import E062, Operator, deliver, endpoints  # noqa: E402

OCTOBER_ONLY = ("nadig_hepg2", "nadig_jurkat", "rpe1", "k562_essential")
# same cell line as the fold's TRUTH -> never used there
FORBIDDEN = {"k562_gwps": {"k562_essential"}}
DOSES = (0.05, 0.10, 0.20, 0.34)
SHUFFLE_SEED = 23323


def deranged(res, rng):
    """Whole-row derangement of a residual: same magnitudes, wrong recipients."""
    n = np.asarray(res.mean_effect).shape[0]
    perm = rng.permutation(n)
    for _ in range(50):
        if not (perm == np.arange(n)).any():
            break
        perm = rng.permutation(n)
    return E062.CD4Prepared(
        mode=res.mode + "|deranged",
        mean_effect=np.asarray(res.mean_effect)[perm].copy(),
        gene_votes=np.asarray(res.gene_votes)[perm].copy(),
        target_condition_count=np.asarray(res.target_condition_count)[perm].copy(),
        condition_details=[])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="*", default=["real", "pess"])
    ap.add_argument("--folds", nargs="*", default=list(P.CANCER_FOLDS))
    ap.add_argument("--sources", nargs="*", default=list(OCTOBER_ONLY))
    ap.add_argument("--full-at", type=float, default=None,
                    help="also compute all six members at this dose")
    ap.add_argument("--out", default=str(P.WORK / "oct233.json"))
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

            def run(arm, extra, g, full=False):
                if (panel, held, arm) in done:
                    return
                if extra is None or g == 0.0:
                    add, per_unit, per_unit_allrows = base_resid, 0.0, 0.0
                else:
                    fx, per_unit = compose218.solve_fraction_carrying(whitened, extra, g)
                    _, per_unit_allrows = plolo.solve_fraction(whitened, extra, g)
                    add = compose218.merged(base_resid, 1.0, extra, fx)
                final, _ = E062.add_post_whiten_residual(whitened, add, 1.0, wcols)
                final[:, fold["panel_positions"]] = 0.0
                m = S.score_arm(final, fold, bctx, gate, ctrl, Ts, LR, pools, real_de,
                                np.random.default_rng(7), None, full=full)
                carrying = int((np.asarray(extra.target_condition_count) > 0).sum()) \
                    if extra is not None else None
                rows.append(dict(panel=panel, fold=held, arm=arm, dose=g,
                                 n_targets=int(final.shape[0]),
                                 carrying_rows=carrying,
                                 per_unit_carrying=per_unit,
                                 per_unit_all_rows=per_unit_allrows,
                                 full=full, seconds=round(time.time() - t0, 1), **m))
                out_path.write_text(json.dumps(rows, indent=2))
                print(f"[{panel}/{held}/{arm}] cos {m['cos_agg']:.6f} pds {m['pds']:.6f}"
                      + (f" carry {carrying} per_unit {per_unit:.4f} "
                         f"(all-rows {per_unit_allrows:.4f})" if extra is not None else ""),
                      flush=True)
                del final

            run("BASE", None, 0.0, full=args.full_at is not None)
            for name in args.sources:
                if name in FORBIDDEN.get(held, set()):
                    print(f"[{panel}/{held}] SKIP {name}: same cell line as the truth",
                          flush=True)
                    continue
                res = resid218.source_residual(name, targets, "gene_norm")
                carry = int((np.asarray(res.target_condition_count) > 0).sum())
                if carry == 0:
                    print(f"[{panel}/{held}] {name}: 0 carrying rows — VACUOUS", flush=True)
                    continue
                shuf = deranged(res, np.random.default_rng(SHUFFLE_SEED))
                for g in DOSES:
                    run(f"{name}@{g:g}", res, g,
                        full=(args.full_at is not None and abs(g - args.full_at) < 1e-9))
                    run(f"{name}_shuf@{g:g}", shuf, g)
            del op, whitened, pre
    out_path.write_text(json.dumps(rows, indent=2))
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
