"""236 Phase 0 + Phase 1 — do the FOUR OCTOBER-ONLY SOURCES COMPOSE?

233 priced `nadig_hepg2`, `nadig_jurkat`, `rpe1`, `k562_essential` INDIVIDUALLY on both of
`DEF-RUNBOOK` §5's synthetic disjoint panels and found all four live (+5.5% to +17.7% relative
aggregate cosine, dose-monotone, clean shuffle control).  It never composed them.  225 measured
today's sources saturating at **1.84x the best single at thirteen sources** -- *"the fifth source
buys 10% of what the first bought"*.  Whether the October four saturate the same way is the
number the October plan turns on, and it is not obvious in either direction: all four are
LOW-COVERAGE sources drawing on overlapping biology (two Nadig lines plus an essential-gene
screen), so they could equally be four names for one asset.

THE INSTRUMENT IS 233's, REUSED VERBATIM.  Same panels, same re-fitted per-panel operators
(axis-asserted), same per-panel CD4 banks, same folds, same residual form (218's un-whitened
post-whitening residual, sized on CARRYING rows -- the all-rows median is 0.0000 for all four
and the shipped rule would deliver identically zero), same renderer and scorer read-back, same
seeds.  Phase 0 re-runs 233's own single-source arms through this file and asserts they
reproduce to 1e-12; if they do not, nothing below is readable.

COMPOSITION FORM.  `compose218.union_rows` -- per-gene mean over whichever members carry a
coordinate, union of carried rows.  That is the H1-union precedent already shipped in
production (222) and 225's own stack form, and it holds the TOTAL dose fixed as `k` grows, so
the k-curve is a composition measurement and not a dose sweep in disguise.  All 15 non-empty
subsets are measured, which yields the greedy order, every one of the 4! orders, and the
leave-one-out cost of each member, from one pass.

LEAK AUDIT.  `k562_essential` is a Replogle K562 screen; it is NEVER used on the fold that
holds `k562_gwps` out as truth (233's rule, kept).  Every stack arm carries a whole-row
DERANGED copy of its own union at the same dose -- the null that holds magnitude fixed.

    uv run python stack236.py --panels pess real
"""
from __future__ import annotations

import argparse
import itertools
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
from plolo import E062, Operator  # noqa: E402

OCTOBER_ONLY = ("nadig_hepg2", "nadig_jurkat", "rpe1", "k562_essential")
FORBIDDEN = {"k562_gwps": {"k562_essential"}}       # same cell line + lab as the fold's truth
DOSES_233 = (0.05, 0.10, 0.20, 0.34)                # 233's grid -- the Phase-0 known-answer rows
DOSES = (0.05, 0.10, 0.20, 0.34, 0.45, 0.55)        # extended: a 4-member union may want more
SHUFFLE_SEED = 23323                                # 233's, so the KA shuffle rows reproduce too
SUBSETS = [s for k in range(1, 5) for s in itertools.combinations(OCTOBER_ONLY, k)]


def tag(names) -> str:
    return "+".join(names)


def deranged(res, rng):
    """233's null: whole-row derangement -- same magnitudes, wrong recipients."""
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
    ap.add_argument("--panels", nargs="*", default=["pess", "real"])
    ap.add_argument("--folds", nargs="*", default=list(P.CANCER_FOLDS))
    ap.add_argument("--out", default=str(WORK / "stack236.json"))
    args = ap.parse_args()

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["panel"], r["fold"], r["arm"]) for r in rows}
    resid218, compose218 = S._resid_module()

    # 233's own numbers, for the Phase-0 assertion
    ka = {(r["panel"], r["fold"], r["arm"]): r
          for r in json.loads((R233 / "work/oct233.json").read_text())}
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
        for held in args.folds:
            t0 = time.time()
            fold = fold233.build(panel, held)
            # FALSIFIER 5: the operator is a function of the panel and fails silently.
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

            allowed = [n for n in OCTOBER_ONLY if n not in FORBIDDEN.get(held, set())]
            src = {n: resid218.source_residual(n, targets, "gene_norm") for n in allowed}
            for n, r in list(src.items()):
                if int((np.asarray(r.target_condition_count) > 0).sum()) == 0:
                    print(f"[{panel}/{held}] {n}: 0 carrying rows — VACUOUS", flush=True)
                    src.pop(n)

            # carrying-row overlap: the mechanism behind whatever the saturation number is
            carry = {n: np.asarray(r.target_condition_count) > 0 for n, r in src.items()}
            overlap = {f"{a}|{b}": int((carry[a] & carry[b]).sum())
                       for a in carry for b in carry}

            def run(arm, extra, g):
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
                                np.random.default_rng(7), None, full=False)
                carrying = int((np.asarray(extra.target_condition_count) > 0).sum()) \
                    if extra is not None else None
                rows.append(dict(panel=panel, fold=held, arm=arm, dose=g,
                                 n_targets=int(final.shape[0]), carrying_rows=carrying,
                                 per_unit_carrying=per_unit,
                                 per_unit_all_rows=per_unit_allrows,
                                 seconds=round(time.time() - t0, 1), **m))
                out_path.write_text(json.dumps(rows, indent=2))
                # PHASE 0: assert against 233's stored value where one exists
                k233 = ka.get((panel, held, arm))
                if k233 is not None:
                    dc = abs(k233["cos_agg"] - m["cos_agg"])
                    dp = abs(k233["pds"] - m["pds"])
                    ok = dc < 1e-12 and dp < 1e-12
                    ka_report.append(dict(panel=panel, fold=held, arm=arm,
                                          d_cos=dc, d_pds=dp, ok=bool(ok),
                                          cos_233=k233["cos_agg"], cos_236=m["cos_agg"]))
                    if not ok:
                        raise SystemExit(
                            f"KA FAILED {panel}/{held}/{arm}: dcos {dc:.3e} dpds {dp:.3e} "
                            f"(233 {k233['cos_agg']:.9f} vs {m['cos_agg']:.9f})")
                    print(f"    KA PASS {arm}: == 233 to {max(dc, dp):.1e}", flush=True)
                print(f"[{panel}/{held}/{arm}] cos {m['cos_agg']:.6f} pds {m['pds']:.6f}"
                      + (f" carry {carrying} per_unit {per_unit:.4f}"
                         f" (all-rows {per_unit_allrows:.4f})" if extra is not None else ""),
                      flush=True)
                del final

            run("BASE", None, 0.0)
            for names in SUBSETS:
                if any(n not in src for n in names):
                    continue
                t = tag(names)
                res = (src[names[0]] if len(names) == 1
                       else compose218.union_rows([src[n] for n in names]))
                shuf = deranged(res, np.random.default_rng(SHUFFLE_SEED))
                doses = DOSES if len(names) > 1 else tuple(sorted(set(DOSES) | set(DOSES_233)))
                for g in doses:
                    run(f"{t}@{g:g}", res, g)
                    # falsifier 1: the derangement holds magnitude fixed
                    run(f"{t}_shuf@{g:g}", shuf, g)
                del res, shuf

            rows.append(dict(panel=panel, fold=held, arm="__overlap__", dose=0.0,
                             n_targets=int(len(targets)), carrying_rows=None,
                             overlap=overlap,
                             carrying={n: int(v.sum()) for n, v in carry.items()},
                             cos_agg=0.0, pds=0.0))
            out_path.write_text(json.dumps(rows, indent=2))
            del op, whitened, pre, src, carry

    (WORK / "ka236.json").write_text(json.dumps(ka_report, indent=2))
    n_ok = sum(1 for r in ka_report if r["ok"])
    print(f"PHASE 0: {n_ok}/{len(ka_report)} known-answer rows reproduce 233 to 1e-12")
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
