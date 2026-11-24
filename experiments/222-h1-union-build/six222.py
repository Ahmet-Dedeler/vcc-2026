"""222 INBOX-1 item 1 — the arm through 216's EXACT `nmae` identity, not through a rate card.

`experiments/216-cd4-hard/six216.py` with its 12-arm ladder replaced by three arms and its
`build_final` replaced by the production builder's merged CD4 + H1-union construction.
Everything else -- 188's exact `nmae` identity, the three pools, `key_precision`, the
delivered-amplitude readout, the a15 reference row -- is 216's, imported, not re-derived.

  a15 f=0.10        the `nmae` reference row (216's base), and this script's KNOWN-ANSWER row:
                    its `nmae wald89` and the incumbent's delta must reproduce
                    `experiments/216-cd4-hard/work/six216.json`.
  cd4 f=0.34        the board incumbent `cd4f034`.
  h1union g=0.28    the arm, built through `scripts/make_submission_pm_cd4.py`'s own
                    `prepare_h1_union` and `solve_carrying_row_fraction`.

Also INBOX-1 item 3: CD4's own carrying-row norm quantile ratio, measured and reported, with
CD4's sizing rule untouched.

    uv run python six222.py
"""
from __future__ import annotations

import gc
import importlib.util
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "experiments/203-reopen-125"))
sys.path.insert(0, str(ROOT / "experiments/216-cd4-hard"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, key_precision, topn, umd  # noqa: E402
from plolo import (E062, FOLDS, OFFICIAL, PANEL, PROD_OPS, Operator,  # noqa: E402
                   build_fold, cd4_for, deliver, solve_fraction)

ALPHA, GAMMA, LAM = 0.982, 0.65, 0.0
NMAE_SPAN = 0.60
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
# (name, cd4 residual f, h1 union g)
ARMS = [
    ("a15 f=0.10 SHIPPED", 0.10, None),
    ("cd4 f=0.34 <- INCUMBENT", 0.34, None),
    ("h1union g=0.28", 0.34, 0.28),
]


def load_builder():
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("six222_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["six222_builder"] = module
    spec.loader.exec_module(module)
    return module


MOD = load_builder()


def build_final(pre, op, fold, cd4, resf, union, g):
    """Production's staging: CD4 by the shipped rule, the union by the carrying-row rule,
    merged into ONE add and ONE global preserve (218 section 7.2)."""
    wh = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)
    frac, _ = solve_fraction(wh, cd4, resf)
    if g is None:
        added, applied = cd4, frac
        fx = 0.0
    else:
        fx, _frame = MOD.solve_carrying_row_fraction(wh, union, g)
        vc = np.asarray(cd4.gene_votes) > 0
        vx = np.asarray(union.gene_votes) > 0
        eff = np.where(vc, np.asarray(cd4.mean_effect), 0.0) * frac
        eff += np.where(vx, np.asarray(union.mean_effect), 0.0) * fx
        added = E062.CD4Prepared(
            mode="merged", mean_effect=eff, gene_votes=(vc | vx).astype(np.uint8),
            target_condition_count=np.maximum(
                np.asarray(cd4.target_condition_count),
                np.asarray(union.target_condition_count)).astype(np.uint8),
            condition_details=[])
        applied = 1.0
    final, _ = E062.add_post_whiten_residual(wh, added, applied, fold["wcols"])
    final[:, fold["panel_positions"]] = 0.0
    return final, float(frac), float(fx)


def delivered_norm(final, fold) -> float:
    """`experiments/216-cd4-hard/form216.py`'s `delivered_norm`, verbatim -- including its
    `fold["scored"]` axis, which is what makes the known-answer row below comparable."""
    vals = []
    for ctx in ("A", "B", "C"):
        b = plolo.challenge_bctrl(ctx)
        expr = np.expm1(b)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        gate = (cpm > plolo.MIN_CPM) & (b > 0)
        dl = deliver(final, b, gate, "minnorm")
        vals.append(float(np.median(np.linalg.norm(dl[:, fold["scored"]], axis=1))))
    return float(np.mean(vals))


def main() -> None:
    acc: dict[tuple[str, str], list] = {}
    prec: dict[str, list] = {}
    norms: dict[str, list] = {}
    pds: dict[str, list] = {}
    cd4_quantiles = []

    for held in FOLDS:
        fold = build_fold(held)
        truth_h = np.asarray(fold["truth"], dtype=np.float64)
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        union, _details = MOD.prepare_h1_union(
            fold["targets"], OFFICIAL, PANEL, max_rss_gib=40.0)
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)

        # INBOX-1 item 3 -- MEASURED ONLY, CD4's sizing rule is not touched.
        added = np.where(np.asarray(cd4.gene_votes) > 0,
                         np.asarray(cd4.mean_effect), 0.0)
        rn = np.linalg.norm(added, axis=1)
        live = rn[rn > 0]
        cd4_quantiles.append(dict(
            fold=held, n_rows=int(rn.size), carrying=int(live.size),
            coverage=float(live.size / rn.size),
            median_all_rows=float(np.median(rn)),
            median_carrying=float(np.median(live)),
            # where the all-rows median sits inside the carrying-row distribution
            all_rows_median_percentile_of_carrying=float(
                100.0 * np.mean(live <= np.median(rn))),
            q14_carrying=float(np.quantile(live, 0.14)),
            overdose_now=float(np.median(live) / np.median(rn)),
            overdose_at_q14=float(np.median(live) / np.quantile(live, 0.14)),
        ))

        bA = plolo.challenge_bctrl("A")
        eA = np.expm1(bA)
        gA = (1e6 * eA / float(eA.sum()) > plolo.MIN_CPM) & (bA > 0)
        sc = fold["scored"]
        cpm = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
        LR = to_lfc(truth_h, fold["b_ctrl"], gA)[:, sc]
        pools = {
            "wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm, 0))[None, :], POOL_N),
            "lfc89": topn(np.abs(LR), POOL_N),
            "h89": topn(np.abs(truth_h[:, sc]), POOL_N),
        }
        Ts = truth_h[:, sc]

        for name, resf, g in ARMS:
            fin, frac, fx = build_final(pre, op, fold, cd4, resf, union, g)
            H = deliver(fin, bA, gA, None)
            LP = to_lfc(H, bA, gA)[:, sc]
            for pk, mmask in pools.items():
                acc.setdefault((name, pk), []).append(umd(LP, LR, mmask))
            prec.setdefault(name, []).append(
                key_precision(np.abs(H[:, sc]), H[:, sc], Ts, n=POOL_N))
            # 216's `six216.py` measures amplitude as the median row norm of the NO-mass-fix
            # render on frame A over the operator columns -- NOT `form216.delivered_norm`,
            # which is the `@adel` solver's target.  Matching it is what makes the
            # known-answer row below comparable.
            norms.setdefault(name, []).append(
                float(np.median(np.linalg.norm(H[:, fold["wcols"]], axis=1))))
            # `pds` on the a23 render, the frame 162 priced
            P = deliver(fin, bA, gA, "minnorm")[:, sc]
            A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
            B = Ts / np.maximum(np.linalg.norm(Ts, axis=1, keepdims=True), 1e-300)
            S = A @ B.T
            own = np.diag(S).copy()
            n = S.shape[0]
            pds.setdefault(name, []).append(
                (held, float(np.mean(((S < own[:, None]).sum(1)
                                      + 0.5 * ((S == own[:, None]).sum(1) - 1))
                                     / (n - 1)))))
            del fin, H, LP, P, A, B, S
            gc.collect()
        print(f"[{held}] done  (cd4 frac {frac:.5f}, union fx {fx:.5f})", flush=True)
        del op, pre, cd4, union, truth_h, fold, LR, Ts, pools
        gc.collect()

    base = ARMS[0][0]
    inc = ARMS[1][0]
    ref = {pk: np.concatenate([x[2] for x in acc[(base, pk)]])
           for pk in ("wald89", "lfc89", "h89")}
    out = {}
    b0 = {}
    for pk in ref:
        U = np.concatenate([x[0] for x in acc[(base, pk)]])
        M = np.concatenate([x[1] for x in acc[(base, pk)]])
        b0[pk] = (float(np.mean(1 + (U - 2 * M) / ref[pk])), U, M)
    p0 = float(np.mean(prec[base]))
    a15n = float(np.mean(norms[base]))

    print("\n" + "=" * 120)
    print("188's EXACT `nmae` identity, ordering precision, delivered amplitude "
          "-- 4 PLOLO folds, challenge frame A")
    print("=" * 120)
    print(f"{'arm':<26}{'dlv/a15':>9}{'nmae wald89':>13}{'d vs a15':>11}"
          f"{'ovr nmae':>11}{'top89 prec':>12}{'d prec':>10}{'nmae lfc89':>12}"
          f"{'nmae h89':>11}{'pds 4f':>9}{'pds 3f':>9}")
    for name, *_ in ARMS:
        vals = {}
        for pk in ref:
            U = np.concatenate([x[0] for x in acc[(name, pk)]])
            M = np.concatenate([x[1] for x in acc[(name, pk)]])
            vals[pk] = (float(np.mean(1 + (U - 2 * M) / ref[pk])), U, M)
        nm = vals["wald89"][0]
        ov = -(nm - b0["wald89"][0]) / NMAE_SPAN / 6.0
        pr = float(np.mean(prec[name]))
        dn = float(np.mean(norms[name])) / a15n
        p4 = float(np.mean([v for _f, v in pds[name]]))
        p3 = float(np.mean([v for f, v in pds[name] if f in CANCER]))
        out[name] = dict(nmae_wald89=nm, d_nmae_vs_a15=nm - b0["wald89"][0],
                         overall_nmae_vs_a15=ov, top89_precision=pr,
                         d_precision_vs_a15=pr - p0, delivered_norm_ratio=dn,
                         nmae_lfc89=vals["lfc89"][0], nmae_h89=vals["h89"][0],
                         pds_4folds=p4, pds_3cancer=p3,
                         pds_per_fold={f: v for f, v in pds[name]})
        print(f"{name:<26}{dn:>9.4f}{nm:>13.5f}{nm-b0['wald89'][0]:>+11.5f}"
              f"{ov:>+11.6f}{pr:>12.5f}{pr-p0:>+10.5f}{vals['lfc89'][0]:>12.5f}"
              f"{vals['h89'][0]:>11.5f}{p4:>9.5f}{p3:>9.5f}")

    ni, pi, di = (out[inc]["nmae_wald89"], out[inc]["top89_precision"],
                  out[inc]["delivered_norm_ratio"])
    print(f"\nsame table, deltas against the board incumbent `{inc}`:")
    print(f"{'arm':<26}{'d nmae':>10}{'ovr nmae':>11}{'d prec':>10}{'d dlv':>9}"
          f"{'d pds 4f':>10}{'d pds 3f':>10}")
    for name, *_ in ARMS:
        o = out[name]
        o["d_nmae_vs_incumbent"] = o["nmae_wald89"] - ni
        o["overall_nmae_vs_incumbent"] = -(o["nmae_wald89"] - ni) / NMAE_SPAN / 6.0
        o["d_precision_vs_incumbent"] = o["top89_precision"] - pi
        o["d_pds_4folds_vs_incumbent"] = o["pds_4folds"] - out[inc]["pds_4folds"]
        o["d_pds_3cancer_vs_incumbent"] = o["pds_3cancer"] - out[inc]["pds_3cancer"]
        print(f"{name:<26}{o['d_nmae_vs_incumbent']:>+10.5f}"
              f"{o['overall_nmae_vs_incumbent']:>+11.6f}"
              f"{o['d_precision_vs_incumbent']:>+10.5f}"
              f"{o['delivered_norm_ratio']-di:>+9.4f}"
              f"{o['d_pds_4folds_vs_incumbent']:>+10.6f}"
              f"{o['d_pds_3cancer_vs_incumbent']:>+10.6f}")

    # KNOWN-ANSWER: 216's own six216.json for the two arms it also ran.
    ka_path = ROOT / "experiments/216-cd4-hard/work/six216.json"
    ka_report = []
    if ka_path.is_file():
        ka = json.loads(ka_path.read_text())
        pairs = [("a15 f=0.10 SHIPPED", "a15 f=0.10 SHIPPED"),
                 ("cd4 f=0.34 <- INCUMBENT", "res f=0.34 <- INCUMBENT")]
        print("\nKNOWN-ANSWER against experiments/216-cd4-hard/work/six216.json:")
        for mine, theirs in pairs:
            if theirs not in ka:
                continue
            for key in ("nmae_wald89", "top89_precision", "delivered_norm_ratio"):
                got, want = out[mine][key], ka[theirs][key]
                ok = abs(got - want) <= 5e-6
                ka_report.append(dict(arm=mine, key=key, got=got, want=want, ok=bool(ok)))
                print(f"  {'PASS' if ok else 'FAIL'}  {mine:26s} {key:22s} "
                      f"got {got:.8f}  216 {want:.8f}  |d| {abs(got-want):.2e}")

    print("\nINBOX-1 item 3 — CD4's own carrying-row norm distribution (MEASURED, not changed):")
    print(f"  {'fold':16s}{'carry':>7}{'cov':>7}{'med all':>10}{'med carry':>11}"
          f"{'pctile':>8}{'now x':>8}{'at q0.14 x':>12}")
    for row in cd4_quantiles:
        print(f"  {row['fold']:16s}{row['carrying']:>7d}{row['coverage']:>7.3f}"
              f"{row['median_all_rows']:>10.4f}{row['median_carrying']:>11.4f}"
              f"{row['all_rows_median_percentile_of_carrying']:>8.1f}"
              f"{row['overdose_now']:>8.4f}{row['overdose_at_q14']:>12.4f}")

    (WORK / "six222.json").write_text(json.dumps(
        dict(arms=out, known_answer=ka_report, cd4_quantiles=cd4_quantiles), indent=2))
    print(f"\nwrote {WORK/'six222.json'}")


if __name__ == "__main__":
    main()
