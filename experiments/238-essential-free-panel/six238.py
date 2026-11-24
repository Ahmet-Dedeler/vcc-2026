"""238 Phase 0.2 + Phase 2 — the four October sources, re-priced, on all six members.

This is `experiments/236-october-stack/six236.py` with three changes and nothing else:

  1. the panel registry is 238's (`bind238`), so `efp` / `efr` -- 199's two draws with the
     common-essential list removed from the eligible pool -- are addressable;
  2. `--ka` asserts every arm that 236 also ran reproduces its stored `cos_agg` / `pds` to
     1e-12.  That is Phase 0.2 of the mandate: if 236's headline does not reproduce through
     this file, nothing below is readable and the run stops;
  3. a source with ZERO carrying rows is recorded as VACUOUS rather than scored.  On an
     essential-free panel that is the expected outcome for `k562_essential`, whose own
     perturbation list IS the exclusion list, and it must be reported as an absence rather
     than as a number that happens to equal BASE.

Everything else -- the residual form (218's un-whitened post-whitening residual, sized on
CARRYING rows), the union composition (`compose218.union_rows`), the whole-row derangement
null, the renderer, the scorer read-back, 188's `e = dM/dU`, the leave-top-3-out falsifier,
and every seed (7 / 23323 / 233233 / 9090+n) -- is 233's and 236's, untouched.

    uv run python six238.py --panels pess real --stacks ... --shuf --ka
    uv run python six238.py --panels efp efr  --stacks ... --doses 0.05 0.1 0.2 0.34 --shuf
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

import bind238  # noqa: E402,F401   MUST be first: rebinds panels233 -> 238's registry
import panels238 as P  # noqa: E402
import score233 as S  # noqa: E402
import fold233  # noqa: E402
import plolo  # noqa: E402
from plolo import E062, Operator, deliver  # noqa: E402

from six236 import leave_top_out, prep_rows  # noqa: E402
from stack236 import FORBIDDEN, SHUFFLE_SEED, deranged  # noqa: E402

R233W = ROOT / "experiments/233-october-retention/work"
# every file holding a (panel, fold, arm) row this experiment can be checked against
KA_FILES = (ROOT / "experiments/236-october-stack/work/six236.json",
            ROOT / "experiments/236-october-stack/work/stack236.json",
            R233W / "score233.json",          # BASE / null / cd4_off / h1_off / resid_off
            R233W / "oct233.json")            # the four sources, 233's own single-source rows
E236 = KA_FILES[0]
E236_STACK = KA_FILES[1]
ESSENTIAL = {ln.strip() for ln in
             (ROOT / "data/donors/k562_essential/perts.txt").read_text().splitlines()
             if ln.strip()}
KA_TOL = 1e-12


def restrict(res, keep_mask):
    """`compose218.restrict_rows` with an explicit row mask: keep those rows, zero the rest.

    Used for the WITHIN-PANEL decomposition.  Comparing `pess` with an essential-free draw
    changes the panel AND the targets at once; restricting the SAME residual on the SAME panel
    to its essential vs its non-essential carrying rows changes only which rows carry it, and
    `solve_fraction_carrying` re-solves the dose on whichever rows survive, so the per-row
    magnitude is matched between the two halves.
    """
    eff = np.asarray(res.mean_effect).copy()
    votes = np.asarray(res.gene_votes).copy()
    counts = np.asarray(res.target_condition_count).copy()
    eff[~keep_mask] = 0.0
    votes[~keep_mask] = 0
    counts[~keep_mask] = 0
    return E062.CD4Prepared(mode=res.mode + "|restricted", mean_effect=eff, gene_votes=votes,
                            target_condition_count=counts, condition_details=[])


def source_perts(name: str) -> set[str]:
    from sources218 import SOURCES
    p = ROOT / SOURCES[name]["path"] / "perts.txt"
    return {ln.strip() for ln in p.read_text().splitlines() if ln.strip()}


def ka_table() -> dict:
    out = {}
    for f in KA_FILES:
        if f.exists():
            for r in json.loads(f.read_text()):
                if r.get("arm") == "__overlap__" or "cos_agg" not in r:
                    continue
                out.setdefault((r["panel"], r["fold"], r["arm"]), r)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="*", default=["efp", "efr"])
    ap.add_argument("--folds", nargs="*", default=list(P.CANCER_FOLDS))
    ap.add_argument("--stacks", nargs="*", default=[],
                    help="'+'-joined source names, no dose (doses come from --doses); "
                         "empty runs BASE and its deranged null only")
    ap.add_argument("--doses", nargs="*", type=float, default=[0.34])
    ap.add_argument("--shuf", action="store_true")
    ap.add_argument("--split-essential", action="store_true",
                    help="also run each stack restricted to its ESSENTIAL carrying rows and "
                         "to its NON-ESSENTIAL ones — the within-panel decomposition")
    ap.add_argument("--ablations", action="store_true",
                    help="also run cd4_off / h1_off / resid_off — 233 §4.1's component "
                         "attribution, and the POSITIVE CONTROL that the instrument still "
                         "resolves a real source on this panel")
    ap.add_argument("--ka", action="store_true",
                    help="assert every arm 236 also ran reproduces to 1e-12")
    ap.add_argument("--leave-top3-at", type=float, default=0.34)
    ap.add_argument("--out", default=str(WORK / "six238.json"))
    args = ap.parse_args()

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["panel"], r["fold"], r["arm"]) for r in rows}
    resid218, compose218 = S._resid_module()
    ka = ka_table()
    ka_report: list[dict] = []

    ctrl_mod = __import__("ctrl228")
    ctrl = {c: ctrl_mod.control_block(c) for c in S.CONTEXTS}
    bctx = {c: plolo.challenge_bctrl(c) for c in S.CONTEXTS}
    gate = {}
    for c in S.CONTEXTS:
        e = np.expm1(bctx[c])
        gate[c] = (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (bctx[c] > 0)

    for panel in args.panels:
        S.bind_panel(panel, resid218)
        panel_set = set(P.read_panel(panel))
        for held in args.folds:
            t0 = time.time()
            fold = fold233.build(panel, held)
            # FALSIFIER 1: the operator is a function of the panel and fails silently.
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
                """FALSIFIER 2 is recorded here: `per_unit_all_rows` is production's own
                median-over-EVERY-row sizing, which returns EXACTLY 0.0 below 50% coverage
                (223's silent zero).  Every dose below is solved on CARRYING rows instead."""
                if extra is None or g == 0.0:
                    add, fx, pu_c, pu_a = base_resid, 0.0, 0.0, 0.0
                else:
                    fx, pu_c = compose218.solve_fraction_carrying(whitened, extra, g)
                    _, pu_a = plolo.solve_fraction(whitened, extra, g)
                    add = compose218.merged(base_resid, 1.0, extra, fx)
                final, _ = E062.add_post_whiten_residual(whitened, add, 1.0, wcols)
                final[:, fold["panel_positions"]] = 0.0
                return final, float(fx), float(pu_c), float(pu_a)

            def score(arm, final, fx, carrying, dose, extra_fields=None):
                m = S.score_arm(final, fold, bctx, gate, ctrl, Ts, LR, pools, real_de,
                                np.random.default_rng(7), [140, 150, 250], full=True)
                H_A = deliver(final, bctx["A"], gate["A"], None)
                LP = to_lfc(H_A, bctx["A"], gate["A"])[:, sc]
                u, mm, d = umd(LP, LR, pools["wald89"])
                m["u_over_d"] = float(np.mean(u / d))
                m["m_over_d"] = float(np.mean(mm / d))
                row = dict(panel=panel, fold=held, arm=arm, dose=dose,
                           n_targets=int(final.shape[0]), carrying_rows=carrying,
                           fraction=fx, seconds=round(time.time() - t0, 1),
                           **(extra_fields or {}), **m)
                rows.append(row)
                out_path.write_text(json.dumps(rows, indent=2))
                k = ka.get((panel, held, arm))
                if k is not None:
                    dc = abs(k["cos_agg"] - m["cos_agg"])
                    dp = abs(k["pds"] - m["pds"])
                    ok = dc < KA_TOL and dp < KA_TOL
                    ka_report.append(dict(panel=panel, fold=held, arm=arm, d_cos=dc,
                                          d_pds=dp, ok=bool(ok), cos_236=k["cos_agg"],
                                          cos_238=m["cos_agg"]))
                    print(f"    KA {'PASS' if ok else '**FAIL**'} {arm}: "
                          f"max|d| {max(dc, dp):.1e} vs 236/233", flush=True)
                    if not ok and args.ka:
                        raise SystemExit(
                            f"KA FAILED {panel}/{held}/{arm}: dcos {dc:.3e} dpds {dp:.3e}")
                print(f"[{panel}/{held}/{arm}] cos {m['cos_agg']:.6f} pds {m['pds']:.6f} "
                      f"nmae {m['nmae']:.5f} fid {m['fid']:.5f} reach {m['reach']:.5f} "
                      f"jac89 {m['jac']['89']:.5f}", flush=True)
                return row

            # ---------------- BASE (kept for leave-top-3-out) and the deranged NULL
            base_final, _, _, _ = build(None, 0.0)
            Pb = {c: prep_rows(deliver(base_final, bctx[c], gate[c], "minnorm")[:, sc], Ts)
                  for c in S.CONTEXTS}
            if (panel, held, "BASE") not in done:
                score("BASE", base_final, 0.0, None, 0.0)
            if (panel, held, "null") not in done:
                rng = np.random.default_rng(S.NULL_SEED)
                n = base_final.shape[0]
                perm = rng.permutation(n)
                for _ in range(50):
                    if not (perm == np.arange(n)).any():
                        break
                    perm = rng.permutation(n)
                score("null", base_final[perm], 0.0, None, 0.0)
            del base_final

            # ---- 233 §4.1's ablations.  `cd4_off` / `h1_off` are the POSITIVE CONTROL: the
            # instrument must still resolve a REAL source on this panel, or a null reading for
            # the October four would be a property of the panel rather than of those sources.
            if args.ablations:
                abl = {
                    "cd4_off": compose218.merged(cd4, 0.0, h1, frac_u),
                    "h1_off": compose218.merged(cd4, frac_cd4, h1, 0.0),
                    "resid_off": compose218.merged(cd4, 0.0, h1, 0.0),
                }
                for nm, add in abl.items():
                    if (panel, held, nm) in done:
                        continue
                    f2, _ = E062.add_post_whiten_residual(whitened, add, 1.0, wcols)
                    f2[:, fold["panel_positions"]] = 0.0
                    score(nm, f2, 0.0, None, 0.0)
                    done.add((panel, held, nm))
                    del f2

            for spec in args.stacks:
                names = tuple(spec.split("+"))
                if any(n in FORBIDDEN.get(held, set()) for n in names):
                    kept = tuple(n for n in names if n not in FORBIDDEN.get(held, set()))
                    print(f"[{panel}/{held}] {spec}: dropping "
                          f"{sorted(set(names) - set(kept))} (same cell line as the truth)",
                          flush=True)
                    names = kept
                    if not names:
                        continue
                srcs = {n: resid218.source_residual(n, targets, "gene_norm") for n in names}
                carry_n = {n: int((np.asarray(r.target_condition_count) > 0).sum())
                           for n, r in srcs.items()}
                vacuous = [n for n, c in carry_n.items() if c == 0]
                if vacuous:
                    for n in vacuous:
                        rec = dict(panel=panel, fold=held, arm=f"{n}@VACUOUS", dose=None,
                                   n_targets=len(targets), carrying_rows=0, vacuous=True,
                                   panel_overlap=len(panel_set & source_perts(n)))
                        if (panel, held, rec["arm"]) not in done:
                            rows.append(rec)
                            done.add((panel, held, rec["arm"]))
                            print(f"[{panel}/{held}] {n}: 0 carrying rows — VACUOUS "
                                  f"(instrument refuses to price it)", flush=True)
                    names = tuple(n for n in names if n not in vacuous)
                    srcs = {n: srcs[n] for n in names}
                    if not names:
                        out_path.write_text(json.dumps(rows, indent=2))
                        continue
                res_all = (srcs[names[0]] if len(names) == 1
                           else compose218.union_rows([srcs[n] for n in names]))
                carrying = int((np.asarray(res_all.target_condition_count) > 0).sum())
                carry_mask = np.asarray(res_all.target_condition_count) > 0
                carry_targets = [targets[i] for i in np.flatnonzero(carry_mask)]
                extra = dict(
                    sources=list(names),
                    carrying_targets_essential=int(sum(t in ESSENTIAL for t in carry_targets)),
                    carrying_targets_essential_frac=(
                        sum(t in ESSENTIAL for t in carry_targets) / carrying)
                    if carrying else None,
                    per_source_carrying=carry_n)
                tag = "+".join(names)
                ess_mask = np.array([t in ESSENTIAL for t in targets], dtype=bool)
                for g in args.doses:
                    plain = f"{tag}@{g:g}"
                    todo = [(plain, res_all, carrying, extra)]
                    if args.shuf:
                        todo.append((f"{tag}_shuf@{g:g}",
                                     deranged(res_all, np.random.default_rng(SHUFFLE_SEED)),
                                     carrying, extra))
                    if args.split_essential and abs(g - args.leave_top3_at) < 1e-9:
                        for lab, mask in (("essonly", carry_mask & ess_mask),
                                          ("noness", carry_mask & ~ess_mask)):
                            n_keep = int(mask.sum())
                            if n_keep == 0:
                                print(f"[{panel}/{held}] {tag}_{lab}@{g:g}: 0 rows — VACUOUS",
                                      flush=True)
                                continue
                            ex = dict(extra, restricted_to=lab,
                                      carrying_targets_essential=(n_keep if lab == "essonly"
                                                                  else 0),
                                      carrying_targets_essential_frac=(
                                          1.0 if lab == "essonly" else 0.0))
                            todo.append((f"{tag}_{lab}@{g:g}",
                                         restrict(res_all, mask), n_keep, ex))
                    for nm, r, n_carry, ex in todo:
                        if (panel, held, nm) in done:
                            continue
                        final, fx, pu_c, pu_a = build(r, g)
                        row = score(nm, final, fx, n_carry, g,
                                    dict(ex, per_unit_carrying=pu_c,
                                         per_unit_all_rows=pu_a))
                        done.add((panel, held, nm))
                        if nm == plain and abs(g - args.leave_top3_at) < 1e-9:
                            lt = {}
                            for c in S.CONTEXTS:
                                Pa = prep_rows(
                                    deliver(final, bctx[c], gate[c], "minnorm")[:, sc], Ts)
                                lt[c] = leave_top_out(Pb[c], Pa)
                                del Pa
                            row["leave_top3"] = {
                                k2: float(np.mean([lt[c][k2] for c in S.CONTEXTS]))
                                for k2 in lt["A"]}
                            rows[-1] = row
                            out_path.write_text(json.dumps(rows, indent=2))
                            print(f"    leave-top-3-out: cos gain "
                                  f"{row['leave_top3']['cos_gain_full']:+.6f} -> "
                                  f"{row['leave_top3']['cos_gain_minus_top3']:+.6f}; "
                                  f"pds gain {row['leave_top3']['pds_gain_full']:+.6f} -> "
                                  f"{row['leave_top3']['pds_gain_minus_top3']:+.6f}",
                                  flush=True)
                        del final
                del srcs, res_all
            del op, whitened, pre, Pb

    out_path.write_text(json.dumps(rows, indent=2))
    if ka_report:
        (WORK / "ka238.json").write_text(json.dumps(ka_report, indent=2))
        n_ok = sum(1 for r in ka_report if r["ok"])
        print(f"\nKNOWN-ANSWER: {n_ok}/{len(ka_report)} arms reproduce 236/233 to {KA_TOL:g}")
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
