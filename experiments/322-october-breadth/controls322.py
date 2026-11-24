"""322 Phase 4 — the mandate's two controls: is the contribution COVERAGE or IDENTITY?

The retention screen (`breadth322.py`) measures what each source is worth on each panel.  It
cannot by itself say whether a source is worth what it is worth BECAUSE it reaches more
targets or because of what it says about them.  Two controls separate that, on 233's panels:

  1. COVERAGE-MATCHED RANDOM ASSIGNMENT.  Hold the source fixed and cut its carrying rows to
     a chosen count at random (3 seeds).  If value tracks coverage, the curve over k is the
     whole story and two different sources at matched k read the same.

  2. COVERAGE HELD FIXED.  Restrict two different sources to the SAME target set (the rows
     both carry on that panel) and compare.  Any difference that survives is biology, not
     reach; any difference that vanishes was reach.

The two live sources with comparable, non-vacuous coverage on ALL THREE panels are
`dld1_lowmoi` (67 / 67 / 73) and `igvf_cardiomyocyte` (87 / 60 / 36).  `kolf2` cannot enter
either control on a disjoint panel: it carries nothing there (`census322.py`).

Everything downstream of the residual is `breadth322.run()`'s path, unchanged.

    uv run python controls322.py --panels abc real pess
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib322 as L  # noqa: E402

import panels233 as P  # noqa: E402
import score233 as S  # noqa: E402
import fold233  # noqa: E402
import plolo  # noqa: E402
from plolo import E062, Operator  # noqa: E402

PAIR = ("dld1_lowmoi", "igvf_cardiomyocyte")
DOSE = 0.34
SEEDS = (322001, 322002, 322003)
LADDER = (20, 40, 60)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="*", default=["abc", "real", "pess"])
    ap.add_argument("--folds", nargs="*", default=list(P.CANCER_FOLDS))
    ap.add_argument("--out", default=str(L.WORK / "controls322.json"))
    args = ap.parse_args()

    ok, msg = L.memory_ok()
    print(f"[mem] {msg} -> {'OK' if ok else 'HOLD'}", flush=True)
    assert ok, "memory gate: hold"

    L.register_extra()
    out_path = Path(args.out)
    rows_out = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["panel"], r["fold"], r["arm"]) for r in rows_out}
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
            cpm_sc = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
            LR = to_lfc(np.asarray(fold["truth"], dtype=np.float64), fold["b_ctrl"],
                        gate["A"])[:, sc]
            pools = {"wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm_sc, 0)), POOL_N),
                     "lfc89": topn(np.abs(LR), POOL_N), "h89": topn(np.abs(Ts), POOL_N)}
            real_de = {n: topn(np.abs(Ts), n) for n in S.N_REALS}

            def run(arm, extra, g, note=None):
                if (panel, held, arm) in done:
                    return
                if extra is None:
                    add, per_unit = base_resid, 0.0
                else:
                    fx, per_unit = compose218.solve_fraction_carrying(whitened, extra, g)
                    add = compose218.merged(base_resid, 1.0, extra, fx)
                final, _ = E062.add_post_whiten_residual(whitened, add, 1.0, wcols)
                final[:, fold["panel_positions"]] = 0.0
                m = S.score_arm(final, fold, bctx, gate, ctrl, Ts, LR, pools, real_de,
                                np.random.default_rng(7), None, full=False)
                carrying = int((np.asarray(extra.target_condition_count) > 0).sum()) \
                    if extra is not None else None
                rows_out.append(dict(panel=panel, fold=held, arm=arm, dose=g,
                                     n_targets=int(final.shape[0]), carrying_rows=carrying,
                                     per_unit_carrying=per_unit, note=note,
                                     seconds=round(time.time() - t0, 1), **m))
                out_path.write_text(json.dumps(rows_out, indent=2))
                print(f"[{panel}/{held}/{arm}] cos {m['cos_agg']:.6f} pds {m['pds']:.6f} "
                      f"carry {carrying}", flush=True)
                del final

            res = {n: resid218.source_residual(n, targets, "gene_norm") for n in PAIR}
            car = {n: L.carrying_rows(res[n]) for n in PAIR}
            common = np.intersect1d(car[PAIR[0]], car[PAIR[1]])
            print(f"[{panel}/{held}] carrying {PAIR[0]} {car[PAIR[0]].size} · "
                  f"{PAIR[1]} {car[PAIR[1]].size} · common {common.size}", flush=True)

            run("CTRL_BASE", None, 0.0)
            if common.size == 0:
                print(f"[{panel}/{held}] no common rows — controls undefined", flush=True)
                del op, whitened, pre
                continue

            # ---- control 2: coverage held fixed (the SAME targets, two biologies) --------
            for n in PAIR:
                run(f"{n}|common{common.size}@{DOSE:g}",
                    L.restrict_residual(res[n], common), DOSE,
                    note="restricted to the rows BOTH sources carry on this panel")

            # ---- control 1: coverage-matched RANDOM restriction (same count, random rows)
            for n in PAIR:
                for s in SEEDS:
                    run(f"{n}|rand{common.size}_s{s}@{DOSE:g}",
                        L.coverage_matched_random(res[n], common.size,
                                                  np.random.default_rng(s)), DOSE,
                        note="random rows of the same source, count matched to `common`")

            # ---- the coverage LADDER on the broader-in-source member ---------------------
            n0 = PAIR[0]
            for k in LADDER:
                if k >= car[n0].size:
                    continue
                for s in SEEDS[:2]:
                    run(f"{n0}|rand{k}_s{s}@{DOSE:g}",
                        L.coverage_matched_random(res[n0], k, np.random.default_rng(s)),
                        DOSE, note="coverage ladder")
            del op, whitened, pre
    out_path.write_text(json.dumps(rows_out, indent=2))
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
