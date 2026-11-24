"""322 Phase 2 — does adding a source to the stack retain its value on a DISJOINT panel?

FORM.  `oct233.py`'s, unchanged: 218's un-whitened post-whitening residual added on top of
the shipped stack (CD4 @ 0.34 + H1 union @ 0.28), at a dose `g` solved on the source's
CARRYING rows, scored on 233's PLOLO folds after the renderer and the scorer's read-back.
Same operators, same CD4 banks, same folds -- all of them 233's files, read-only.

WHY THIS FORM.  233 §7 item 3's argument is about the RESIDUAL COMPONENTS of the stack:
"the component whose contribution collapses is the one whose coverage collapses (CD4,
24.7 % -> 10.7 %), and the one that is coverage-robust by construction (the H1 union, sized
on carrying rows) does not move."  CD4 and the H1 union are both post-whitening residual
sources sized on rows.  So a candidate "breadth" purchase is measured in the same form.
264's alternative -- `kolf2` as a fourth BASE DONOR -- is handled in `basedonor322.py`.

THE UNIT.  The resampled unit is the PANEL DRAW: 233 built two (`real`, `pess`) and the
control (`abc`).  n = 2 disjoint draws.  Every "retention" number here is a ratio of two
measured contributions and its uncertainty is bounded by that draw count -- three folds are
NOT three draws (they share the panel).

LEAK AUDIT.  A source is never used on the fold that holds its own cell line out as truth
(`FORBIDDEN`).  A whole-row-deranged copy of each source is run at the same dose as the null.

    uv run python breadth322.py --panels abc real pess
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
from oct233 import deranged  # noqa: E402  (233's whole-row derangement, reused verbatim)

# The candidates.  Each row is (name, why it is here).
CANDIDATES = (
    "kolf2",              # 273/300 on TODAY's panel, 273 perturbations TOTAL  <- the mandate's
    "dld1_lowmoi",        # 67/300 today, 4,715 perturbations TOTAL            <- the contrast
    "telohaec",           # 65/300 today, 2,167 total
    "igvf_cardiomyocyte",  # 87/300 today, 1,888 total
)
# same cell line as the fold's TRUTH -> never used there.  233's rule; no candidate here is
# K562, HCT116 or HEK293T, so the set is empty, but the guard stays so it is auditable.
FORBIDDEN: dict[str, set[str]] = {"k562_gwps": {"k562_essential"}}
DOSES = (0.20, 0.34)
FULL_AT = 0.34
SHUFFLE_SEED = 23323


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="*", default=["abc", "real", "pess"])
    ap.add_argument("--folds", nargs="*", default=list(P.CANCER_FOLDS))
    ap.add_argument("--sources", nargs="*", default=list(CANDIDATES))
    ap.add_argument("--out", default=str(L.WORK / "breadth322.json"))
    ap.add_argument("--register", nargs="*", default=[],
                    help="NAME=PATH_TEMPLATE, e.g. kolf2_refit=work/kolf2-refit/{panel}; "
                         "the template is formatted with the panel key so a source that must "
                         "be re-extracted per panel can enter the same screen")
    args = ap.parse_args()

    ok, msg = L.memory_ok()
    print(f"[mem] {msg} -> {'OK' if ok else 'HOLD'}", flush=True)
    assert ok, "memory gate: hold"

    L.register_extra()
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

    import sources218
    for panel in args.panels:
        S.bind_panel(panel, resid218)
        for spec in args.register:
            nm, tmpl = spec.split("=", 1)
            p = (L.HERE / tmpl.format(panel=panel)).resolve()
            assert (p / "perts.txt").exists(), f"{nm}: no bank at {p}"
            entry = dict(path=str(p.relative_to(L.ROOT)),
                         role=f"per-panel re-extraction ({panel})")
            # `resid218` did `from sources218 import SOURCES`, so it holds a REFERENCE to
            # the dict object.  Mutate both bindings; rebinding `sources218.SOURCES` after
            # `resid218` is imported would leave the importer looking at the old dict.
            sources218.SOURCES[nm] = entry
            resid218.SOURCES[nm] = entry
            resid218._BANK_CACHE.pop(nm, None)
            print(f"[register] {nm} -> {p.relative_to(L.ROOT)} "
                  f"({sum(1 for _ in (p / 'perts.txt').open())} perturbations)", flush=True)
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

            def run(arm, extra, g, full=False, note=None):
                if (panel, held, arm) in done:
                    return
                if extra is None or g == 0.0:
                    add, per_unit, per_unit_allrows = base_resid, 0.0, 0.0
                else:
                    fx, per_unit = compose218.solve_fraction_carrying(whitened, extra, g)
                    # 233 falsifier 4 / 223's silent zero.  Below 50 % coverage the
                    # production ALL-ROWS median solve returns exactly 0.0 -- the residual
                    # is not small, it is identically absent.  252 turned that into a
                    # refusal, so the zero now arrives as an exception; record it as the
                    # 0.0 it is (that IS the measurement 233 §5.1 reports).
                    try:
                        _, per_unit_allrows = plolo.solve_fraction(whitened, extra, g)
                    except plolo.SilentZeroError:
                        per_unit_allrows = 0.0
                    add = compose218.merged(base_resid, 1.0, extra, fx)
                final, _ = E062.add_post_whiten_residual(whitened, add, 1.0, wcols)
                final[:, fold["panel_positions"]] = 0.0
                m = S.score_arm(final, fold, bctx, gate, ctrl, Ts, LR, pools, real_de,
                                np.random.default_rng(7), None, full=full)
                carrying = int((np.asarray(extra.target_condition_count) > 0).sum()) \
                    if extra is not None else None
                rows.append(dict(panel=panel, fold=held, arm=arm, dose=g,
                                 n_targets=int(final.shape[0]), carrying_rows=carrying,
                                 per_unit_carrying=per_unit,
                                 per_unit_all_rows=per_unit_allrows, full=full, note=note,
                                 seconds=round(time.time() - t0, 1), **m))
                out_path.write_text(json.dumps(rows, indent=2))
                print(f"[{panel}/{held}/{arm}] cos {m['cos_agg']:.6f} pds {m['pds']:.6f}"
                      + (f" carry {carrying} per_unit {per_unit:.4f} "
                         f"(all-rows {per_unit_allrows:.4f})" if extra is not None else ""),
                      flush=True)
                del final

            run("BASE", None, 0.0, full=True)
            for name in args.sources:
                if name in FORBIDDEN.get(held, set()):
                    print(f"[{panel}/{held}] SKIP {name}: same cell line as the truth",
                          flush=True)
                    continue
                res = resid218.source_residual(name, targets, "gene_norm")
                carry = int((np.asarray(res.target_condition_count) > 0).sum())
                if carry == 0:
                    # 223's SILENT ZERO, made loud.  This is not a skipped row: it is the
                    # measurement.  A source with no carrying rows contributes EXACTLY
                    # nothing while still being named in `VCC_UNION_SOURCES`.
                    rows.append(dict(panel=panel, fold=held, arm=f"{name}@VACUOUS", dose=None,
                                     n_targets=len(targets), carrying_rows=0,
                                     per_unit_carrying=0.0, per_unit_all_rows=0.0,
                                     full=False, note="0 carrying rows — contributes exactly 0",
                                     seconds=round(time.time() - t0, 1)))
                    out_path.write_text(json.dumps(rows, indent=2))
                    print(f"[{panel}/{held}] {name}: 0 carrying rows — VACUOUS", flush=True)
                    continue
                shuf = deranged(res, np.random.default_rng(SHUFFLE_SEED))
                for g in DOSES:
                    run(f"{name}@{g:g}", res, g, full=abs(g - FULL_AT) < 1e-9)
                run(f"{name}_shuf@{FULL_AT:g}", shuf, FULL_AT)
            del op, whitened, pre
    out_path.write_text(json.dumps(rows, indent=2))
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
