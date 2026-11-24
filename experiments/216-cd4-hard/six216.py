"""216 Phase 3 — the six members at every candidate operating point.

Route, and why each member takes the route it does (rule 12: price all six or name the
mechanism that pins the ones you cannot):

  `pds`   PLOLO, 4 folds, a23 render.  Priceable (162: sign 6/6, LOO sigma 0.001494 raw).
  `mse`   096/197's identity through `X` only, in `price216.py`.  **211's standing rule
          applies: `claim` is at 100.0000% of its cap in the one context that reads it, so
          NO arm here may be re-priced on amplitude grounds** -- which is exactly why the
          library form's alpha-doubling cannot be sold as an `mse` gain.
  `nmae`  188's EXACT identity on the four PLOLO folds with the `|lfc|*sqrt(cpm)` pool at the
          reference's own BH depth (89 genes/row) -- the construction 188 validated on
          h1-test.  Both alternative pools are carried, because 188 showed the construction
          can invert the sign.
  `reach` / `fid` / `jac`   **UNMEASURED on h1-test, with the mechanism**: `vacuity216.py`
          re-measures that the CD4 residual is IDENTICALLY zero on h1-test (median row norm
          0.0 over 85 rows, 81 of them empty) so h1-test cannot see this arm at all.  The
          ORDERING half of `reach`/`fid` is readable on PLOLO and is reported; the yield
          halves carry `n_real` and are unreadable on every harness we hold (rule 9).

    uv run python six216.py
"""
from __future__ import annotations

import gc
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
from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   deliver, solve_fraction)
from form216 import build_final, delivered_norm, solve_alpha_delivered, votes_for  # noqa: E402
from sweep216 import ALPHA, GAMMA, LAM  # noqa: E402

NMAE_SPAN = 0.60

# (name, library weight, residual f, alpha policy)   policy: "fixed" | "adel"
ARMS = [
    ("a15 f=0.10 SHIPPED", 0.0, 0.10, "fixed"),
    ("res f=0.20", 0.0, 0.20, "fixed"),
    ("res f=0.28", 0.0, 0.28, "fixed"),
    ("res f=0.34 <- INCUMBENT", 0.0, 0.34, "fixed"),
    ("res f=0.40", 0.0, 0.40, "fixed"),
    ("res f=0.55", 0.0, 0.55, "fixed"),
    ("res f=0.80", 0.0, 0.80, "fixed"),
    ("res f=1.00 <- cos argmax", 0.0, 1.00, "fixed"),
    ("lib w=1 (1/k) @a0", 1.0, 0.0, "fixed"),
    ("lib w=1 (1/k) @adel", 1.0, 0.0, "adel"),
    ("lib w=0.5 @adel", 0.5, 0.0, "adel"),
    ("lib w=0.5 + res f=0.34 @adel", 0.5, 0.34, "adel"),
]


def main() -> None:
    acc: dict[tuple[str, str], list] = {}
    prec: dict[str, list] = {}
    norms: dict[str, list] = {}
    alphas: dict[str, list] = {}

    for held in FOLDS:
        fold = build_fold(held)
        truth_h = np.asarray(fold["truth"], dtype=np.float64)
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        cd4e = np.asarray(cd4.mean_effect, dtype=np.float64)
        carry = np.asarray(cd4.target_condition_count) == 3
        votes = votes_for(held, fold["targets"])
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        w = fold["wcols"]

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

        inc_final, _ = build_final(pre, ALPHA, op, fold, cd4, 0.10)
        target_norm = delivered_norm(inc_final, fold)
        del inc_final
        gc.collect()

        for name, lw, rf, policy in ARMS:
            sel = carry & (votes > 0)
            pre4 = pre.copy()
            if lw > 0:
                pre4[sel] = ((pre[sel] * votes[sel, None] + lw * cd4e[sel])
                             / (votes[sel, None] + lw))
            a = ALPHA
            if policy == "adel":
                a, _ = solve_alpha_delivered(pre4, op, fold, cd4, rf, target_norm)
            fin, _ = build_final(pre4, a, op, fold, cd4, rf)
            H = deliver(fin, bA, gA, None)
            LP = to_lfc(H, bA, gA)[:, sc]
            for pk, mmask in pools.items():
                acc.setdefault((name, pk), []).append(umd(LP, LR, mmask))
            prec.setdefault(name, []).append(
                key_precision(np.abs(H[:, sc]), H[:, sc], Ts, n=POOL_N))
            norms.setdefault(name, []).append(
                float(np.median(np.linalg.norm(H[:, w], axis=1))))
            alphas.setdefault(name, []).append(a)
            del fin, H, LP, pre4
            gc.collect()
        print(f"[{held}] done", flush=True)
        del op, pre, cd4, truth_h, fold, LR, Ts, pools
        gc.collect()

    base = ARMS[0][0]
    inc = "res f=0.34 <- INCUMBENT"
    out = {}
    ref = {pk: np.concatenate([x[2] for x in acc[(base, pk)]])
           for pk in ("wald89", "lfc89", "h89")}
    b0 = {}
    for pk in ref:
        U = np.concatenate([x[0] for x in acc[(base, pk)]])
        M = np.concatenate([x[1] for x in acc[(base, pk)]])
        b0[pk] = (float(np.mean(1 + (U - 2 * M) / ref[pk])), U, M)
    p0 = float(np.mean(prec[base]))
    a15n = float(np.mean(norms[base]))

    print(f"\n{'='*132}")
    print("`nmae` (188 exact identity), ordering precision, delivered amplitude — "
          "4 PLOLO folds, challenge frame A")
    print("=" * 132)
    print(f"{'arm':<30}{'alpha':>7}{'dlv/a15':>9}{'nmae wald89':>13}{'d':>10}{'e':>8}"
          f"{'ovr nmae':>11}{'top89 prec':>12}{'d prec':>10}{'nmae lfc89':>12}"
          f"{'nmae h89':>11}")
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
        al = float(np.mean(alphas[name]))
        dn = float(np.mean(norms[name])) / a15n
        out[name] = dict(nmae_wald89=nm, d_nmae_vs_a15=nm - b0["wald89"][0], e=e,
                         overall_nmae_vs_a15=ov, top89_precision=pr,
                         d_precision_vs_a15=pr - p0, alpha=al, delivered_norm_ratio=dn,
                         nmae_lfc89=vals["lfc89"][0], nmae_h89=vals["h89"][0])
        print(f"{name:<30}{al:>7.3f}{dn:>9.4f}{nm:>13.5f}{nm-b0['wald89'][0]:>+10.5f}"
              f"{e:>8.3f}{ov:>+11.6f}{pr:>12.5f}{pr-p0:>+10.5f}"
              f"{vals['lfc89'][0]:>12.5f}{vals['h89'][0]:>11.5f}")

    # everything again, relative to the CURRENT incumbent rather than to a15
    print(f"\nsame table, deltas taken against the incumbent `{inc}`:")
    print(f"{'arm':<30}{'d nmae':>10}{'ovr nmae':>11}{'d prec':>10}{'d dlv':>9}")
    ni, pi, di = out[inc]["nmae_wald89"], out[inc]["top89_precision"], \
        out[inc]["delivered_norm_ratio"]
    for name, *_ in ARMS:
        o = out[name]
        o["d_nmae_vs_incumbent"] = o["nmae_wald89"] - ni
        o["overall_nmae_vs_incumbent"] = -(o["nmae_wald89"] - ni) / NMAE_SPAN / 6.0
        o["d_precision_vs_incumbent"] = o["top89_precision"] - pi
        print(f"{name:<30}{o['d_nmae_vs_incumbent']:>+10.5f}"
              f"{o['overall_nmae_vs_incumbent']:>+11.6f}"
              f"{o['d_precision_vs_incumbent']:>+10.5f}"
              f"{o['delivered_norm_ratio']-di:>+9.4f}")
    (WORK / "six216.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {WORK/'six216.json'}")


if __name__ == "__main__":
    main()
