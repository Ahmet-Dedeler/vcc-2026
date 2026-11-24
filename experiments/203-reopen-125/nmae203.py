"""203 — `nmae` at each arm's OWN re-solved amplitude.

`nmae` is the one member that is genuinely magnitude-sensitive (187: a15's alpha is the
joint `mse`+`nmae` optimum and the surface is a sharp rho cliff), so pricing a 4.7-6.8x
shrunk estimator at a15's alpha measures the shrink.  `phase23.py` priced every arm at
a15's alpha; this prices the same arms at the amplitude that matches a15's delivered
median row norm, which is what the mandate's "re-solve the amplitude" step asks for.

KNOWN-ANSWER ROW: `a15` at alpha* == `a15` at alpha (alpha* is 0.982 by construction), so
its `nmae*` must reproduce `phase23.json`'s wald89 a15 value exactly.
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
from phase23 import POOL_N, topn, umd  # noqa: E402
from plolo import FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver, solve_fraction  # noqa: E402

OUT = HERE / "work"
ALPHA, GAMMA, LAM, CD4F = 0.982, 0.65, 0.0, 0.10
SHRINK, NMAE_SPAN = 6.7801, 0.60


def mrn(M, w):
    return float(np.median(np.linalg.norm(np.asarray(M)[:, w], axis=1))) or 1.0


def main():
    acc, meta = {}, {}
    for held in FOLDS:
        fold = build_fold(held, limit_gib=7.0)
        truth_h = np.asarray(fold["truth"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                      if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
        z = np.load(OUT / f"sources-{held}.npz", allow_pickle=True)
        panel = est203.prepare(np.asarray(z["effects"], dtype=np.float64),
                               available=z["available"], support=z["support"])
        targets = [str(t) for t in z["targets"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        w = fold["wcols"]
        n0 = mrn(pre, w)

        eq, _, _ = est203.equal_mean_agreement(panel)
        arms = {"a15": (pre, ALPHA),
                "a15@alpha* (KA, == a15)": (pre, ALPHA * n0 / mrn(pre, w)),
                "eq122": (eq, ALPHA),
                "eq122@alpha*": (eq, ALPHA * n0 / mrn(eq, w)),
                f"shrink{SHRINK:.2f}": (pre / SHRINK, ALPHA),
                f"shrink{SHRINK:.2f}@alpha*": (pre / SHRINK, ALPHA * SHRINK)}
        for nm in ("eb_broad_prior", "eb_conservative"):
            e, _, _ = est203.nested_posterior(panel, targets, est203.SETTINGS[nm])
            arms[f"post125:{nm}"] = (e, ALPHA)
            arms[f"post125:{nm}@alpha*"] = (e, ALPHA * n0 / mrn(e, w))

        bA = plolo.challenge_bctrl("A")
        eA = np.expm1(bA)
        gA = (1e6 * eA / float(eA.sum()) > plolo.MIN_CPM) & (bA > 0)
        sc = fold["scored"]
        cpm = 1e6 * np.expm1(fold["b_ctrl"][sc]) / float(np.expm1(fold["b_ctrl"]).sum())
        LR = to_lfc(truth_h, fold["b_ctrl"], gA)[:, sc]
        pool = topn(np.abs(LR) * np.sqrt(np.maximum(cpm, 0))[None, :], POOL_N)

        for name, (eff, a) in arms.items():
            p = np.asarray(eff, dtype=np.float64) * a
            wh = op.apply(p, w, GAMMA, LAM)
            frac, _ = solve_fraction(wh, cd4, CD4F)
            fin, _ = plolo.E062.add_post_whiten_residual(wh, cd4, frac, w)
            fin[:, fold["panel_positions"]] = 0.0
            H = deliver(fin, bA, gA, None)
            acc.setdefault(name, []).append(umd(to_lfc(H, bA, gA)[:, sc], LR, pool))
            meta.setdefault(name, []).append(dict(alpha=float(a),
                                                  delivered_norm=mrn(H, w)))
        print(f"[{held}] done", flush=True)
        del op

    b = acc["a15"]
    D0 = np.concatenate([x[2] for x in b])
    U0 = np.concatenate([x[0] for x in b])
    M0 = np.concatenate([x[1] for x in b])
    nm0 = float(np.mean(1 + (U0 - 2 * M0) / D0))
    print(f"\n{'='*104}\n`nmae` (188 identity, lfc space, wald89 significance pool), "
          f"each arm at a15's alpha AND at its own alpha*\n{'='*104}")
    print(f"{'arm':<32}{'alpha':>8}{'dlv norm/a15':>14}{'nmae*':>10}{'d nmae*':>11}"
          f"{'e':>8}{'d overall':>12}")
    out = {}
    a15n = float(np.mean([m["delivered_norm"] for m in meta["a15"]]))
    for name in acc:
        v = acc[name]
        U = np.concatenate([x[0] for x in v])
        M = np.concatenate([x[1] for x in v])
        nm = float(np.mean(1 + (U - 2 * M) / D0))
        du = float(np.mean((U - U0) / D0))
        dm = float(np.mean((M - M0) / D0))
        e = dm / du if abs(du) > 1e-12 else float("nan")
        ov = -(nm - nm0) / NMAE_SPAN / 6.0
        al = float(np.mean([m["alpha"] for m in meta[name]]))
        dn = float(np.mean([m["delivered_norm"] for m in meta[name]])) / a15n
        out[name] = dict(nmae=nm, dnmae=nm - nm0, e=e, overall=ov, alpha=al,
                         delivered_norm_ratio=dn)
        print(f"{name:<32}{al:>8.3f}{dn:>14.4f}{nm:>10.5f}{nm-nm0:>+11.5f}{e:>8.3f}"
              f"{ov:>+12.6f}")
    (OUT / "nmae203.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {OUT/'nmae203.json'}")


if __name__ == "__main__":
    main()
