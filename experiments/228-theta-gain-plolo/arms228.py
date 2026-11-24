"""228 Phase 2 — the gain as a per-gene multiplier on the PRODUCTION effect matrix, priced
on 162's PLOLO instrument after `profile_from_effect` + the scorer's read-back.

WIRING (preregistered)
----------------------
    pre --(x ALPHA*a)--> Operator(gamma=0.65, lam=0) --> whitened
    whitened --(x m, then median-row-norm preserve over wcols)--> whitened'
    whitened' + merged(CD4 @ f=0.34, H1 union @ g=0.28) --> final, panel coords zeroed
    final --> deliver (a23 minnorm and a15 raw renders, challenge frames A/B/C)

    m = 1 + t * (g / median(g) - 1)      t = 0 is the shipped stack, EXACTLY.

The preserve rescale is production's own convention (`Operator.apply`,
`add_post_whiten_residual`), so at fixed `a` the arm is amplitude-neutral by construction and
the `VCC_ALPHA` interaction is a SEPARATE axis, swept here as `a`.

KNOWN-ANSWER ROWS -- the run aborts if either fails:
  KA-1  t=0, CD4 only          == 223 `focus223` BASE:cd4      per fold, to 1e-12
  KA-2  t=0, CD4 + H1 union    == 223 `focus223` BASE:base222  per fold, to 1e-12

    uv run python arms228.py --folds orion_hct116
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
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   solve_fraction)

sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
from resid218 import ALPHA, GAMMA, LAM, score_effects, source_residual  # noqa: E402
from compose218 import merged, solve_fraction_carrying, union_rows  # noqa: E402

CD4_F, UNION_G = 0.34, 0.28
UNION = ("gse295214_h1", "h1_2025")

# 223 `work/focus223.json`, the shipped stack's own numbers -- the known-answer rows
KA = {
    ("k562_gwps", "cd4"): (0.03635247771108033, 0.7958811888512827),
    ("k562_gwps", "base222"): (0.03943565390646537, 0.8053968077826482),
    ("orion_hct116", "cd4"): (0.02968920878381447, 0.790140396929563),
    ("orion_hct116", "base222"): (0.032080093772601896, 0.7922125442469077),
    ("orion_hek293t", "cd4"): (0.024558952922866917, 0.7605504090145728),
    ("orion_hek293t", "base222"): (0.026156872210626163, 0.7636228552854588),
    ("kolf2", "cd4"): (0.009016327626527735, 0.6145272211448684),
    ("kolf2", "base222"): (0.012520242539233312, 0.6361551030668678),
}

ARM_DOSES = {
    "gain": (0.05, 0.10, 0.15, 0.25, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0),
    "gain_reg": (0.05, 0.10, 0.15, 0.25, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0),
    "gain_fullcov": (0.15, 0.25, 0.5, 1.0, 1.5),
    "shuf": (0.15, 0.25, 0.5, 1.0, 1.5),
    "mm": (0.15, 0.25, 0.5, 1.0, 1.5),
    "derange": (0.15, 0.25, 0.5, 1.0, 1.5),
    "basal12": (0.15, 0.25, 0.5, 1.0, 1.5),
    "basal12_wrong": (0.15, 0.25, 0.5, 1.0, 1.5),
}
ALPHA_SCALES = (0.85, 0.925, 1.0, 1.075, 1.15)
ALPHA_ARMS = ("none", "gain", "gain_reg")


def multiplier(g: np.ndarray, live: np.ndarray, t: float) -> np.ndarray:
    med = float(np.median(g[live]))
    assert abs(med) > 1e-12, "degenerate gain median"
    return 1.0 + t * (g / med - 1.0)


def apply_gain(whitened: np.ndarray, m: np.ndarray, wcols) -> tuple[np.ndarray, float]:
    """Multiply, then restore the base's median row norm over the operator columns."""
    out = whitened * m
    before = float(np.median(np.linalg.norm(whitened[:, wcols], axis=1))) or 1.0
    after = float(np.median(np.linalg.norm(out[:, wcols], axis=1))) or 1.0
    out *= before / after
    return out, before / after


def run_fold(held: str, args) -> list[dict]:
    t0 = time.time()
    fold = build_fold(held)
    wcols = np.asarray(fold["wcols"])
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    pre = np.asarray(fold["pre"], dtype=np.float64)
    targets = fold["targets"]
    cd4 = cd4_for(targets)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    op_path = (PROD_OPS / "whitener_union_pruned.npz"
               if FOLDS[held]["op_kind"] == "production" else Path(FOLDS[held]["op"]))
    op = Operator(op_path)

    z = np.load(WORK / f"gain-{held}.npz", allow_pickle=True)
    live = np.asarray(z["live"], dtype=bool)
    gains = {k: np.asarray(z[k], dtype=np.float64) for k in
             ("gain", "gain_reg", "gain_fullcov", "shuf", "mm", "derange",
              "basal12", "basal12_wrong")}

    rows: list[dict] = []
    out_path = WORK / f"arms228-{held}.json"
    if out_path.exists():
        rows = json.loads(out_path.read_text())
    done = {(r["arm"], r["a"]) for r in rows}

    def price(name: float, a: float, m):
        whitened = op.apply(pre * (ALPHA * a), wcols, GAMMA, LAM)
        scale = 1.0
        if m is not None:
            whitened, scale = apply_gain(whitened, m, wcols)
        frac_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
        frac_u, _ = solve_fraction_carrying(whitened, h1, UNION_G)
        res = merged(cd4, frac_cd4, h1, frac_u)
        final, _ = E062.add_post_whiten_residual(whitened, res, 1.0, wcols)
        final[:, fold["panel_positions"]] = 0.0
        e = score_effects(final, fold, truth)
        e = {k: v for k, v in e.items() if not k.startswith("_")}
        # CD4-only variant is only needed for the known-answer row
        return e, dict(frac_cd4=float(frac_cd4), frac_union=float(frac_u),
                       gain_preserve_scale=float(scale),
                       median_row_norm_whitened=float(
                           np.median(np.linalg.norm(whitened[:, wcols], axis=1))))

    # ---------------- KA rows: t = 0 must reproduce the shipped stack exactly
    if ("BASE", 1.0) not in done:
        whitened = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
        frac_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
        frac_u, _ = solve_fraction_carrying(whitened, h1, UNION_G)
        for tag, res in (("cd4", merged(cd4, frac_cd4, None, 0.0)),
                         ("base222", merged(cd4, frac_cd4, h1, frac_u))):
            final, _ = E062.add_post_whiten_residual(whitened, res, 1.0, wcols)
            final[:, fold["panel_positions"]] = 0.0
            e = score_effects(final, fold, truth)
            exp_cos, exp_pds = KA[(held, tag)]
            dc, dp = abs(e["cos"] - exp_cos), abs(e["pds"] - exp_pds)
            print(f"[{held}] KA {tag:8s} cos {e['cos']:.15f} (d {dc:.2e})  "
                  f"pds {e['pds']:.15f} (d {dp:.2e})", flush=True)
            assert dc < 1e-12 and dp < 1e-12, f"KA FAILED for {held}/{tag}"
            if tag == "base222":
                rows.append(dict(fold=held, arm="BASE", variant="none", t=0.0, a=1.0,
                                 **{k: v for k, v in e.items() if not k.startswith("_")},
                                 frac_cd4=float(frac_cd4), frac_union=float(frac_u),
                                 gain_preserve_scale=1.0))
        del whitened
        out_path.write_text(json.dumps(rows, indent=2))

    # ---------------- the dose grid at production alpha
    for name, doses in ARM_DOSES.items():
        for t in doses:
            arm = f"{name}@{t:g}"
            if (arm, 1.0) in done:
                continue
            e, aud = price(arm, 1.0, multiplier(gains[name], live, t))
            rows.append(dict(fold=held, arm=arm, variant=name, t=float(t), a=1.0,
                             **e, **aud))
            out_path.write_text(json.dumps(rows, indent=2))
            print(f"[{held}/{arm}] cos_agg {e['cos_agg']:.6f} pds {e['pds']:.6f} "
                  f"({time.time() - t0:.0f}s)", flush=True)

    # ---------------- the VCC_ALPHA interaction
    for name in ALPHA_ARMS:
        for a in ALPHA_SCALES:
            arm = f"{name}@1" if name != "none" else "none@0"
            if (arm, a) in done or (a == 1.0 and (arm, 1.0) in
                                    {(r["arm"], r["a"]) for r in rows}):
                continue
            m = None if name == "none" else multiplier(gains[name], live, 1.0)
            e, aud = price(arm, a, m)
            rows.append(dict(fold=held, arm=arm, variant=name,
                             t=0.0 if name == "none" else 1.0, a=float(a), **e, **aud))
            out_path.write_text(json.dumps(rows, indent=2))
            print(f"[{held}/{arm} a={a:g}] cos_agg {e['cos_agg']:.6f} "
                  f"pds {e['pds']:.6f} ({time.time() - t0:.0f}s)", flush=True)

    out_path.write_text(json.dumps(rows, indent=2))
    print(f"[{held}] {len(rows)} arms in {time.time() - t0:.0f}s -> {out_path}", flush=True)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    args = ap.parse_args()
    for held in args.folds:
        run_fold(held, args)


if __name__ == "__main__":
    main()
