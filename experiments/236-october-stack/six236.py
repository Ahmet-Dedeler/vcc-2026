"""236 Phase 2 — the SIX MEMBERS of the best October stack, plus 188's `e = dM/dU`.

WHY THIS EXISTS, in two prior results.  228's arm survived read-back on cosine AND on `pds`
and then died at **-0.00702** on the six-member composite, because the cosine it bought was
bought with MAGNITUDES: `e = dM/dU = 0.072` against a break-even of 0.5, so 93% of the
magnitude it moved onto the reference's DE pool was unmatched and `nmae` charged for it.
229's `ARM6` read `+0.00268` on a rate card and measured `-0.00219`.  AGENTS.md 8.1d is
therefore explicit: **before pricing any direction arm, answer whether it moves SIGNS or
MAGNITUDES** -- and the statistic that answers it is `e`, computed through 188's identity on
the `wald89` pool that 188 validated on h1-test.

A source ADDITION should read as SIGNS (`h1union`'s precedent: all six members rose,
`nmae` +0.00588); a rescale reads as MAGNITUDES.  This file measures which one the October
stack is, on the same folds, panels and read-back as Phase 1.

FALSIFIER 4 (mandate) — LEAVE-TOP-3-OUT RETENTION, on `cos_agg` and `pds` SEPARATELY.  These
sources carry only 25-48 rows of a 150-300 target fold, so a gain concentrated in a handful of
targets would be an artefact of three lucky rows rather than an asset.  Each target's
influence is its own leave-one-out effect on the GAIN (arm minus BASE, recomputed with that
row dropped from both, and -- for `pds` -- dropped from the decoy set too, which is the
honest form for a retrieval statistic).  The three most influential are removed and the
retained fraction of the gain reported.

    uv run python six236.py --stacks "rpe1@0.34" "nadig_hepg2+nadig_jurkat+rpe1@0.34"
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
WORK.mkdir(parents=True, exist_ok=True)
R233 = ROOT / "experiments/233-october-retention"
sys.path.insert(0, str(R233))

import panels233 as P  # noqa: E402
import score233 as S  # noqa: E402
import fold233  # noqa: E402
import plolo  # noqa: E402
from plolo import E062, Operator, deliver  # noqa: E402

from stack236 import FORBIDDEN, SHUFFLE_SEED, deranged  # noqa: E402

TOP_OUT = 3


def prep_rows(Pm, Tm):
    """Everything both endpoints need, precomputed ONCE so a leave-one-out is O(n^2), not a
    fresh (n x 17.9k) matmul per row."""
    x = (Pm * Tm).sum(1)
    a = (Pm * Pm).sum(1)
    d = (Tm * Tm).sum(1)
    A = Pm / np.maximum(np.linalg.norm(Pm, axis=1, keepdims=True), 1e-300)
    B = Tm / np.maximum(np.linalg.norm(Tm, axis=1, keepdims=True), 1e-300)
    return dict(x=x, a=a, d=d, S=A @ B.T)


def cos_from(pr, keep):
    """Identical algebra to `plolo.endpoints`' `sqrt(1 - sse)`, restricted to rows `keep`."""
    x, a, d = pr["x"][keep].sum(), pr["a"][keep].sum(), pr["d"][keep].sum()
    return float(x / np.sqrt(a * d)) if a > 0 and d > 0 else 0.0


def pds_from(pr, keep):
    """`pds` restricted to `keep` -- which drops those rows from the DECOY SET too, the honest
    form for a retrieval statistic."""
    s = pr["S"][np.ix_(keep, keep)]
    n = s.shape[0]
    own = np.diag(s).copy()
    beaten = (s < own[:, None]).sum(1)
    ties = (s == own[:, None]).sum(1) - 1
    return float(np.mean((beaten + 0.5 * ties) / (n - 1)))


def leave_top_out(prb, pra, k=TOP_OUT):
    """Per-row influence on the GAIN, then the gain with the `k` most influential dropped."""
    n = prb["x"].shape[0]
    idx = np.arange(n)
    g_cos = cos_from(pra, idx) - cos_from(prb, idx)
    g_pds = pds_from(pra, idx) - pds_from(prb, idx)
    inf_cos = np.empty(n)
    inf_pds = np.empty(n)
    for i in range(n):
        m = np.delete(idx, i)
        inf_cos[i] = g_cos - (cos_from(pra, m) - cos_from(prb, m))
        inf_pds[i] = g_pds - (pds_from(pra, m) - pds_from(prb, m))
    out = {}
    for name, inf, gain, fn in (("cos", inf_cos, g_cos, cos_from),
                                ("pds", inf_pds, g_pds, pds_from)):
        keep = np.setdiff1d(idx, np.argsort(-inf)[:k])
        left = fn(pra, keep) - fn(prb, keep)
        out[f"{name}_gain_full"] = float(gain)
        out[f"{name}_gain_minus_top{k}"] = float(left)
        out[f"{name}_retained"] = float(left / gain) if gain != 0 else float("nan")
        out[f"{name}_top{k}_share"] = float(inf[np.argsort(-inf)[:k]].sum() / gain) \
            if gain != 0 else float("nan")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="*", default=["pess", "real"])
    ap.add_argument("--folds", nargs="*", default=list(P.CANCER_FOLDS))
    ap.add_argument("--stacks", nargs="*", required=True,
                    help="e.g. 'rpe1@0.34' 'nadig_hepg2+rpe1@0.34'")
    ap.add_argument("--shuf", action="store_true", help="also run each stack's derangement")
    ap.add_argument("--out", default=str(WORK / "six236.json"))
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
            from phase23 import POOL_N, topn, umd
            cpm_sc = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(
                np.expm1(fold["b_ctrl"]).sum())
            LR = to_lfc(np.asarray(fold["truth"], dtype=np.float64),
                        fold["b_ctrl"], gate["A"])[:, sc]
            pools = {"wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm_sc, 0)), POOL_N),
                     "lfc89": topn(np.abs(LR), POOL_N), "h89": topn(np.abs(Ts), POOL_N)}
            real_de = {n: topn(np.abs(Ts), n) for n in S.N_REALS}

            def build(extra, g):
                if extra is None or g == 0.0:
                    add, fx = base_resid, 0.0
                else:
                    fx, _ = compose218.solve_fraction_carrying(whitened, extra, g)
                    add = compose218.merged(base_resid, 1.0, extra, fx)
                final, _ = E062.add_post_whiten_residual(whitened, add, 1.0, wcols)
                final[:, fold["panel_positions"]] = 0.0
                return final, float(fx)

            def score(arm, final, fx, carrying, dose):
                m = S.score_arm(final, fold, bctx, gate, ctrl, Ts, LR, pools, real_de,
                                np.random.default_rng(7), [150, 250], full=True)
                # 188's identity: U (predicted magnitude on the reference's DE pool) and
                # M (its sign-matched part).  `e = dM/dU` decides SIGNS vs MAGNITUDES.
                H_A = deliver(final, bctx["A"], gate["A"], None)
                LP = to_lfc(H_A, bctx["A"], gate["A"])[:, sc]
                u, mm, d = umd(LP, LR, pools["wald89"])
                m["u_over_d"] = float(np.mean(u / d))
                m["m_over_d"] = float(np.mean(mm / d))
                row = dict(panel=panel, fold=held, arm=arm, dose=dose,
                           n_targets=int(final.shape[0]), carrying_rows=carrying,
                           fraction=fx, seconds=round(time.time() - t0, 1), **m)
                rows.append(row)
                out_path.write_text(json.dumps(rows, indent=2))
                print(f"[{panel}/{held}/{arm}] cos {m['cos_agg']:.6f} pds {m['pds']:.6f} "
                      f"nmae {m['nmae']:.5f} fid {m['fid']:.5f} reach {m['reach']:.5f} "
                      f"jac89 {m['jac']['89']:.5f} U/D {m['u_over_d']:.5f} "
                      f"M/D {m['m_over_d']:.5f}", flush=True)
                return row

            # ---------------- BASE, kept in memory for the leave-top-3-out comparison
            base_final, _ = build(None, 0.0)
            Pb = {c: prep_rows(deliver(base_final, bctx[c], gate[c], "minnorm")[:, sc], Ts)
                  for c in S.CONTEXTS}
            if (panel, held, "BASE") not in done:
                score("BASE", base_final, 0.0, None, 0.0)
            del base_final

            for spec in args.stacks:
                names_s, _, dose_s = spec.partition("@")
                names = tuple(names_s.split("+"))
                g = float(dose_s)
                if any(n in FORBIDDEN.get(held, set()) for n in names):
                    kept = tuple(n for n in names if n not in FORBIDDEN.get(held, set()))
                    print(f"[{panel}/{held}] {spec}: dropping "
                          f"{set(names) - set(kept)} (same cell line as the truth)", flush=True)
                    names = kept
                    if not names:
                        continue
                srcs = [resid218.source_residual(n, targets, "gene_norm") for n in names]
                res = srcs[0] if len(srcs) == 1 else compose218.union_rows(srcs)
                carrying = int((np.asarray(res.target_condition_count) > 0).sum())
                arm = f"{'+'.join(names)}@{g:g}"
                todo = [(arm, res)]
                if args.shuf:
                    todo.append((f"{'+'.join(names)}_shuf@{g:g}",
                                 deranged(res, np.random.default_rng(SHUFFLE_SEED))))
                for nm, r in todo:
                    if (panel, held, nm) in done:
                        continue
                    final, fx = build(r, g)
                    row = score(nm, final, fx, carrying, g)
                    if not nm.endswith(f"_shuf@{g:g}"):
                        lt = {}
                        for c in S.CONTEXTS:
                            Pa = prep_rows(
                                deliver(final, bctx[c], gate[c], "minnorm")[:, sc], Ts)
                            lt[c] = leave_top_out(Pb[c], Pa)
                            del Pa
                        row["leave_top3"] = {
                            k: float(np.mean([lt[c][k] for c in S.CONTEXTS]))
                            for k in lt["A"]}
                        rows[-1] = row
                        out_path.write_text(json.dumps(rows, indent=2))
                        print(f"    leave-top-{TOP_OUT}-out: cos gain "
                              f"{row['leave_top3']['cos_gain_full']:+.6f} -> "
                              f"{row['leave_top3']['cos_gain_minus_top3']:+.6f} "
                              f"({100*row['leave_top3']['cos_retained']:.1f}% retained); "
                              f"pds gain {row['leave_top3']['pds_gain_full']:+.6f} -> "
                              f"{row['leave_top3']['pds_gain_minus_top3']:+.6f} "
                              f"({100*row['leave_top3']['pds_retained']:.1f}%)", flush=True)
                    del final
                del srcs, res
            del op, whitened, pre, Pb
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
