"""207 Phase 3 — the exact `nmae` identity and the ordering half of `reach`/`fid`,
at each candidate operating point AND at its own re-solved alpha.

WHY NOT h1-test, which the mandate names.  MEASURED, not assumed: the CD4 banks on disk
(`data/external/cd4t/cloud-compact/cd4_compact`) are PANEL-RESTRICTED -- 283/287/284
perturbations, all inside the 2026 300-gene panel.  h1-test's 85 targets intersect the panel in
**5 genes**, and all three CD4 conditions cover only **4** of them.  A CD4-dose arm is therefore
a no-op on 81 of h1-test's 85 rows, and h1-test cannot price this arm at all.  `nmae` is priced
here instead by 188's identity on the PLOLO folds with the pool construction that 188 VALIDATED
on h1-test (`|lfc|*sqrt(cpm)` at the reference's BH depth, 89 genes/row); the two alternative
pools are reported as the sensitivity 188 requires, because the construction can invert the sign.

`nmae`  = mean_p ( 1 + (U_p - 2 M_p) / D_p )   over the pooled coordinates (188/169/203):
        U = sum |pred lfc|, M = sum of sign-matched min(|pred|,|truth|), D = sum |truth lfc|.
        D is held at the incumbent's, so the arms are compared on the same denominator.
`reach`/`fid`  ORDERING halves only: |truth|-weighted sign precision of the top-89 coordinates
        under the production key |H|.  The yield / `n_real` halves are unreadable on any local
        harness (rule 9) and are reported as such.

    uv run python six207.py
"""
from __future__ import annotations

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

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, key_precision, topn, umd  # noqa: E402
from plolo import FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver, solve_fraction  # noqa: E402
from library207 import votes_for  # noqa: E402
from sweep207 import ALPHA, GAMMA, LAM  # noqa: E402

NMAE_SPAN = 0.60          # 169/188/203's conservative end of AGENTS.md §4's 0.57-0.63

# (name, library weight, residual f, re-solve alpha?)
ARMS = [
    ("a15  SHIPPED", 0.0, 0.10, False),
    ("res f=0.14", 0.0, 0.14, False),
    ("res f=0.20", 0.0, 0.20, False),
    ("res f=0.24", 0.0, 0.24, False),
    ("res f=0.28", 0.0, 0.28, False),
    ("res f=0.28 @a*", 0.0, 0.28, True),
    ("res f=0.34", 0.0, 0.34, False),
    ("res f=0.34 @a*", 0.0, 0.34, True),
    ("res f=0.40", 0.0, 0.40, False),
    ("res f=0.40 @a*", 0.0, 0.40, True),
    ("lib w=0.5", 0.5, 0.0, True),
    ("lib w=0.5 @alpha0.982", 0.5, 0.0, False),
    ("lib w=0.5 + res f=0.10", 0.5, 0.10, True),
    ("lib w=0.5 + res f=0.20", 0.5, 0.20, True),
    ("lib w=0.375 + res f=0.20", 0.375, 0.20, True),
    ("lib w=0.75 + res f=0.20", 0.75, 0.20, True),
    ("lib w=0.5 + res f=0.28", 0.5, 0.28, True),
]


def main() -> None:
    acc: dict[tuple[str, str], list] = {}
    prec: dict[str, list] = {}
    norms: dict[str, list] = {}
    for held in FOLDS:
        fold = build_fold(held, limit_gib=7.0)
        truth_h = np.asarray(fold["truth"], dtype=np.float64)
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        cd4e = np.asarray(cd4.mean_effect, dtype=np.float64)
        carry = np.asarray(cd4.target_condition_count) == 3
        votes = votes_for(held, fold["targets"])
        op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                      if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
        w = fold["wcols"]
        n0 = float(np.median(np.linalg.norm(pre[:, w], axis=1)))

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

        for name, lw, rf, resolve in ARMS:
            sel = carry & (votes > 0)
            pre4 = pre.copy()
            if lw > 0:
                pre4[sel] = ((pre[sel] * votes[sel, None] + lw * cd4e[sel])
                             / (votes[sel, None] + lw))
            n1 = float(np.median(np.linalg.norm(pre4[:, w], axis=1)))
            a = ALPHA * n0 / n1 if (resolve and n1) else ALPHA
            wh = op.apply(pre4 * a, w, GAMMA, LAM)
            if rf > 0:
                frac, _ = solve_fraction(wh, cd4, rf)
                fin, _ = plolo.E062.add_post_whiten_residual(wh, cd4, frac, w)
            else:
                fin = wh
            fin[:, fold["panel_positions"]] = 0.0
            H = deliver(fin, bA, gA, None)
            LP = to_lfc(H, bA, gA)[:, sc]
            for pk, m in pools.items():
                acc.setdefault((name, pk), []).append(umd(LP, LR, m))
            prec.setdefault(name, []).append(
                key_precision(np.abs(H[:, sc]), H[:, sc], Ts, n=POOL_N))
            norms.setdefault(name, []).append(
                (float(np.median(np.linalg.norm(H[:, w], axis=1))), a))
        print(f"[{held}] done", flush=True)
        del op

    base = ARMS[0][0]
    out = {}
    print(f"\n{'='*128}\n`nmae` (188 identity), ordering precision, and delivered amplitude, "
          f"four PLOLO folds, challenge frame A\n{'='*128}")
    print(f"{'arm':<26}{'alpha':>7}{'dlv/a15':>9}{'nmae wald89':>13}{'Δ':>10}"
          f"{'e':>8}{'ovr nmae':>10}{'top89 prec':>12}{'Δprec':>10}"
          f"{'nmae lfc89':>12}{'nmae h89':>11}")
    ref = {pk: np.concatenate([x[2] for x in acc[(base, pk)]]) for pk in ("wald89", "lfc89", "h89")}
    b0 = {}
    for pk in ref:
        U = np.concatenate([x[0] for x in acc[(base, pk)]])
        M = np.concatenate([x[1] for x in acc[(base, pk)]])
        b0[pk] = (float(np.mean(1 + (U - 2 * M) / ref[pk])), U, M)
    p0 = float(np.mean(prec[base]))
    a15n = float(np.mean([x[0] for x in norms[base]]))
    for name, *_ in ARMS:
        vals = {}
        for pk in ref:
            U = np.concatenate([x[0] for x in acc[(name, pk)]])
            M = np.concatenate([x[1] for x in acc[(name, pk)]])
            vals[pk] = (float(np.mean(1 + (U - 2 * M) / ref[pk])), U, M)
        nm, U, M = vals["wald89"]
        du = float(np.mean((U - b0["wald89"][1]) / ref["wald89"]))
        dm = float(np.mean((M - b0["wald89"][2]) / ref["wald89"]))
        e = dm / du if abs(du) > 1e-12 else float("nan")
        ov = -(nm - b0["wald89"][0]) / NMAE_SPAN / 6.0
        pr = float(np.mean(prec[name]))
        al = float(np.mean([x[1] for x in norms[name]]))
        dn = float(np.mean([x[0] for x in norms[name]])) / a15n
        out[name] = dict(nmae_wald89=nm, d_nmae=nm - b0["wald89"][0], e=e, overall_nmae=ov,
                         top89_precision=pr, d_precision=pr - p0, alpha=al,
                         delivered_norm_ratio=dn,
                         nmae_lfc89=vals["lfc89"][0], nmae_h89=vals["h89"][0])
        print(f"{name:<26}{al:>7.3f}{dn:>9.4f}{nm:>13.5f}{nm-b0['wald89'][0]:>+10.5f}"
              f"{e:>8.3f}{ov:>+10.6f}{pr:>12.5f}{pr-p0:>+10.5f}"
              f"{vals['lfc89'][0]:>12.5f}{vals['h89'][0]:>11.5f}")
    (WORK / "six207.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {WORK/'six207.json'}")


if __name__ == "__main__":
    main()
