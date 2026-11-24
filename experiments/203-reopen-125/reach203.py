"""203 — the ORDERING half of `reach` / `fid` / `jac`, per arm, on its OWN profile.

Rule 12 forbids shipping an arm with unmeasured members that can move.  `reach`, `fid`
and `jac` cannot be priced in board units on any harness we hold -- 162 measures `reach`
at 2/6 correct signs (worse than chance) and `fid`/`jac` at none, and rule 9 explains why:
all three carry `n_real` in a denominator or a union and our harnesses sit in the wrong
regime.  What CAN be measured, and is what those members read, is the ORDERING: given the
reference's significance-selected pool, how well does the arm's own `|H|` ranking recover
it, and how well does it sign it.

  reach*  recall of the reference's wald89 pool inside the arm's own top-89 by |H|
  fid*    |truth|-weighted sign precision inside the arm's own top-89 by |H|
  jac*    Jaccard of the arm's top-89 with the reference's wald89 pool

Reported as ordering proxies on the PLOLO folds, NOT as board deltas.  The `n_real` regime
error (rule 9) has already cost four board predictions and these numbers do not escape it.
"""
from __future__ import annotations

import json
import os
import sys
import warnings
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

warnings.filterwarnings("ignore", message="All-NaN slice encountered")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "src"))

import est203  # noqa: E402
import plolo  # noqa: E402
from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, topn  # noqa: E402
from plolo import FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver, solve_fraction  # noqa: E402

OUT = HERE / "work"
ALPHA, GAMMA, LAM, CD4F, SHRINK = 0.982, 0.65, 0.0, 0.10, 6.7801


def main():
    acc = {}
    for held in FOLDS:
        fold = build_fold(held, limit_gib=7.0)
        truth = np.asarray(fold["truth"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                      if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
        z = np.load(OUT / f"sources-{held}.npz", allow_pickle=True)
        panel = est203.prepare(np.asarray(z["effects"], dtype=np.float64),
                               available=z["available"], support=z["support"])
        targets = [str(t) for t in z["targets"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        eq, _, _ = est203.equal_mean_agreement(panel)
        arms = {"a15": pre, "eq122": eq, f"shrink{SHRINK:.2f}": pre / SHRINK}
        for nm in ("eb_broad_prior", "eb_conservative"):
            arms[f"post125:{nm}"] = est203.nested_posterior(
                panel, targets, est203.SETTINGS[nm])[0]

        bA = plolo.challenge_bctrl("A")
        eA = np.expm1(bA)
        gA = (1e6 * eA / float(eA.sum()) > plolo.MIN_CPM) & (bA > 0)
        sc = fold["scored"]
        cpm = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
        LR = to_lfc(truth, fold["b_ctrl"], gA)[:, sc]
        ref = topn(np.abs(LR) * np.sqrt(np.maximum(cpm, 0))[None, :], POOL_N)
        Ts = truth[:, sc]

        for name, eff in arms.items():
            p = np.asarray(eff, dtype=np.float64) * ALPHA
            wh = op.apply(p, fold["wcols"], GAMMA, LAM)
            frac, _ = solve_fraction(wh, cd4, CD4F)
            fin, _ = plolo.E062.add_post_whiten_residual(wh, cd4, frac, fold["wcols"])
            fin[:, fold["panel_positions"]] = 0.0
            H = deliver(fin, bA, gA, None)[:, sc]
            call = topn(np.abs(H), POOL_N)
            rec, fid, jac = [], [], []
            for i in range(H.shape[0]):
                r = np.flatnonzero(ref[i])
                c = np.flatnonzero(call[i])
                inter = np.intersect1d(r, c).size
                rec.append(inter / max(r.size, 1))
                jac.append(inter / max(np.union1d(r, c).size, 1))
                w = np.abs(Ts[i, c])
                ok = np.sign(H[i, c]) == np.sign(Ts[i, c])
                if w.sum() > 0:
                    fid.append(float(np.where(ok, w, 0.0).sum() / w.sum()))
            acc.setdefault(name, {"reach": [], "fid": [], "jac": []})
            acc[name]["reach"].append(float(np.mean(rec)))
            acc[name]["fid"].append(float(np.mean(fid)))
            acc[name]["jac"].append(float(np.mean(jac)))
        print(f"[{held}] done", flush=True)
        del op

    print(f"\n{'='*88}\nORDERING PROXIES for `reach`/`fid`/`jac` -- PLOLO 4-fold mean, "
          f"pool = wald{POOL_N}\n{'='*88}")
    print(f"{'arm':<30}{'reach*':>10}{'d':>10}{'fid*':>10}{'d':>10}{'jac*':>10}{'d':>10}")
    base = {k: float(np.mean(v)) for k, v in acc["a15"].items()}
    out = {}
    for name, v in acc.items():
        m = {k: float(np.mean(x)) for k, x in v.items()}
        out[name] = dict(m, d_reach=m["reach"] - base["reach"],
                         d_fid=m["fid"] - base["fid"], d_jac=m["jac"] - base["jac"])
        print(f"{name:<30}{m['reach']:>10.5f}{m['reach']-base['reach']:>+10.5f}"
              f"{m['fid']:>10.5f}{m['fid']-base['fid']:>+10.5f}"
              f"{m['jac']:>10.5f}{m['jac']-base['jac']:>+10.5f}")
    print("\nNOT board deltas.  Rule 9: `n_real` is 89 here and ~2,000 on the challenge; "
          "only the\nordering-driven part of these transfers, and 162 scores `reach` at "
          "2/6 signs locally.")
    (OUT / "reach203.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {OUT/'reach203.json'}")


if __name__ == "__main__":
    main()
