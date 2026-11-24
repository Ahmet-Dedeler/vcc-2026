"""228 — leave-top-3-out retention, on cosine and on `pds` SEPARATELY.

`h1union` realised ~100% of its cosine prediction and only 70% of its `pds` one, and 222
measured why: its Δ`pds` was carried by three targets per fold (leave-top-3-out retained
14.6%) while its cosine gain was distributed (72.7%).  `docs/NEXT-BUILD-SPEC.md` §4 makes
this test mandatory for any new arm.

Also here: the paired-target bootstrap on the per-target cosine, which is the interval the
preregistration requires.

    uv run python retain228.py --arm gain_reg@0.25 [--adel]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
from plolo import (FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   deliver)

sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
from resid218 import source_residual  # noqa: E402
from compose218 import union_rows  # noqa: E402

from arms228 import UNION, multiplier  # noqa: E402
from six228 import build_final  # noqa: E402

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
BOOT = 4000
SEED = 20260904


def per_target(final, fold, truth, bctx, gate):
    """Per-target x, a, d (for the energy-weighted aggregate) and per-target pds,
    averaged over the three challenge frames, on the a23 (minnorm) render."""
    sc = fold["scored"]
    X = np.zeros(final.shape[0]); A = np.zeros_like(X); D = np.zeros_like(X)
    P = np.zeros_like(X)
    for c in ("A", "B", "C"):
        dl = deliver(final, bctx[c], gate[c], "minnorm")[:, sc]
        X += (dl * truth).sum(1); A += (dl * dl).sum(1); D += (truth * truth).sum(1)
        u = dl / np.maximum(np.linalg.norm(dl, axis=1, keepdims=True), 1e-300)
        v = truth / np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-300)
        S = u @ v.T
        n = S.shape[0]
        own = np.diag(S).copy()
        beaten = (S < own[:, None]).sum(1)
        ties = (S == own[:, None]).sum(1) - 1
        P += (beaten + 0.5 * ties) / (n - 1)
    return X / 3.0, A / 3.0, D / 3.0, P / 3.0


def agg(x, a, d, keep=None):
    s = slice(None) if keep is None else keep
    return float(x[s].sum() / np.sqrt(max(a[s].sum(), 1e-300) * max(d[s].sum(), 1e-300)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="gain_reg@0.25")
    ap.add_argument("--adel", action="store_true",
                    help="use the amplitude-matched alpha from adel228.json")
    args = ap.parse_args()
    rng = np.random.default_rng(SEED)
    adel = {(r["fold"], r["arm"]): r["a"] for r in
            json.loads((WORK / "adel228.json").read_text())} if args.adel else {}

    bctx = {c: plolo.challenge_bctrl(c) for c in ("A", "B", "C")}
    gate = {}
    for c in ("A", "B", "C"):
        e = np.expm1(bctx[c])
        gate[c] = (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (bctx[c] > 0)

    out = {}
    for held in FOLDS:
        fold = build_fold(held)
        wcols = np.asarray(fold["wcols"])
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        h1 = union_rows([source_residual(n, fold["targets"], "gene_norm") for n in UNION])
        op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                      if FOLDS[held]["op_kind"] == "production" else Path(FOLDS[held]["op"]))
        z = np.load(WORK / f"gain-{held}.npz", allow_pickle=True)
        live = np.asarray(z["live"], dtype=bool)
        name, _, tt = args.arm.partition("@")
        m = multiplier(np.asarray(z[name], dtype=np.float64), live, float(tt))
        a = adel.get((held, args.arm), 1.0)

        fb = build_final(fold, op, pre, cd4, h1, wcols, None, 1.0)
        xb, ab, db, pb = per_target(fb, fold, truth, bctx, gate)
        del fb
        fa = build_final(fold, op, pre, cd4, h1, wcols, m, a)
        xa, aa, da, pa = per_target(fa, fold, truth, bctx, gate)
        del fa, op

        n = len(xb)
        full_b, full_a = agg(xb, ab, db), agg(xa, aa, da)
        d_full = full_a - full_b
        # per-target contribution to the aggregate gain, and to the pds gain
        cb = xb / np.sqrt(np.maximum(ab * db, 1e-300))
        ca = xa / np.sqrt(np.maximum(aa * da, 1e-300))
        dcos = ca - cb
        dpds = pa - pb
        top_c = np.argsort(-dcos)[:3]
        top_p = np.argsort(-dpds)[:3]
        keep_c = np.ones(n, bool); keep_c[top_c] = False
        keep_p = np.ones(n, bool); keep_p[top_p] = False
        ret_cos = (agg(xa, aa, da, keep_c) - agg(xb, ab, db, keep_c)) / d_full
        ret_pds = float(dpds[keep_p].mean() / dpds.mean()) if dpds.mean() != 0 else np.nan
        share_p = float(dpds[top_p].sum() / (n * dpds.mean())) if dpds.mean() != 0 else np.nan

        idx = rng.integers(0, n, size=(BOOT, n))
        bs = np.array([agg(xa[i], aa[i], da[i]) - agg(xb[i], ab[i], db[i]) for i in idx])
        bsp = dpds[idx].mean(1)

        out[held] = dict(
            n_targets=int(n), a=float(a),
            cos_agg_base=full_b, cos_agg_arm=full_a, d_cos_agg=float(d_full),
            rel_cos_agg_pct=100.0 * d_full / abs(full_b),
            d_pds=float(dpds.mean()),
            retention_cos_leave_top3=float(ret_cos),
            retention_pds_leave_top3=float(ret_pds),
            top3_share_of_d_pds=share_p,
            boot_cos_ci=[float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))],
            boot_cos_p_gt0=float((bs > 0).mean()),
            boot_pds_ci=[float(np.quantile(bsp, 0.025)), float(np.quantile(bsp, 0.975))],
            boot_pds_p_gt0=float((bsp > 0).mean()),
            frac_targets_cos_up=float((dcos > 0).mean()),
            frac_targets_pds_up=float((dpds > 0).mean()))
        r = out[held]
        print(f"[{held}] a {a:.4f}  d_cos_agg {d_full:+.6f} ({r['rel_cos_agg_pct']:+.3f}%)  "
              f"retention_cos {ret_cos:.3f}   d_pds {r['d_pds']:+.6f}  "
              f"retention_pds {ret_pds:.3f}  top3 share {share_p:.3f}\n"
              f"          boot cos 95% [{r['boot_cos_ci'][0]:+.6f}, "
              f"{r['boot_cos_ci'][1]:+.6f}] P>0 {r['boot_cos_p_gt0']:.4f}   "
              f"boot pds 95% [{r['boot_pds_ci'][0]:+.6f}, {r['boot_pds_ci'][1]:+.6f}] "
              f"P>0 {r['boot_pds_p_gt0']:.4f}   targets up cos "
              f"{r['frac_targets_cos_up']:.3f} pds {r['frac_targets_pds_up']:.3f}",
              flush=True)

    cm = {k: float(np.mean([out[f][k] for f in CANCER])) for k in
          ("d_cos_agg", "rel_cos_agg_pct", "d_pds", "retention_cos_leave_top3",
           "retention_pds_leave_top3", "top3_share_of_d_pds")}
    out["CANCER_MEAN"] = cm
    print(f"\nCANCER MEAN  rel cos_agg {cm['rel_cos_agg_pct']:+.3f}%  "
          f"d_pds {cm['d_pds']:+.6f}  retention cos {cm['retention_cos_leave_top3']:.3f}  "
          f"retention pds {cm['retention_pds_leave_top3']:.3f}")
    tag = "adel" if args.adel else "a1"
    p = WORK / f"retain228-{args.arm.replace('@', '_')}-{tag}.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"-> {p}")


if __name__ == "__main__":
    main()
