"""228 — 208's rider rule: is the gain's cosine already supplied by the shipped stack?

The headline is measured MARGINAL over `base222` (three donors + PM50 + CD4 `f = 0.34` +
H1 union `g = 0.28`, board `0.2118310`).  208's rule says a rider must be re-priced against
the CURRENT incumbent, never the one it was measured on -- and the diagnostic that goes with
it is whether the marginal SHRINKS as the stack grows.  If the gain's contribution is
unchanged from the bare whitened consensus to the full stack, it is orthogonal to everything
already in it; if it collapses, the stack already supplies it.

    uv run python stack228.py
"""
from __future__ import annotations

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
from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold,  # noqa: E402
                   cd4_for, solve_fraction)

sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
from resid218 import ALPHA, GAMMA, LAM, score_effects, source_residual  # noqa: E402
from compose218 import merged, solve_fraction_carrying, union_rows  # noqa: E402

from arms228 import CD4_F, UNION, UNION_G, apply_gain, multiplier  # noqa: E402

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
DOSE = 0.25


def main() -> None:
    out = {}
    for held in CANCER:
        fold = build_fold(held)
        wcols = np.asarray(fold["wcols"])
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        h1 = union_rows([source_residual(n, fold["targets"], "gene_norm") for n in UNION])
        op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                      if FOLDS[held]["op_kind"] == "production" else Path(FOLDS[held]["op"]))
        z = np.load(WORK / f"gain-{held}.npz", allow_pickle=True)
        m = multiplier(np.asarray(z["gain_reg"], dtype=np.float64),
                       np.asarray(z["live"], bool), DOSE)

        for gain_on in (False, True):
            wh = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
            if gain_on:
                wh, _ = apply_gain(wh, m, wcols)
            fc, _ = solve_fraction(wh, cd4, CD4_F)
            fu, _ = solve_fraction_carrying(wh, h1, UNION_G)
            stacks = {"bare": None,
                      "cd4": merged(cd4, fc, None, 0.0),
                      "base222": merged(cd4, fc, h1, fu)}
            for sname, res in stacks.items():
                f = (wh.copy() if res is None
                     else E062.add_post_whiten_residual(wh, res, 1.0, wcols)[0])
                f[:, fold["panel_positions"]] = 0.0
                e = score_effects(f, fold, truth)
                out[f"{held}|{sname}|{'gain' if gain_on else 'none'}"] = dict(
                    cos_agg=e["cos_agg"], pds=e["pds"])
                del f
            del wh
        del op
        print(f"[{held}] done", flush=True)

    print(f"\n{'MARGINAL of gain_reg@0.25 as the stack grows (a=1, 3 cancer folds)':^78}")
    print(f"{'stack':<12}{'cos_agg base':>14}{'cos_agg +gain':>15}{'rel %':>9}"
          f"{'d_pds':>10}")
    summary = {}
    for sname in ("bare", "cd4", "base222"):
        b = float(np.mean([out[f"{f}|{sname}|none"]["cos_agg"] for f in CANCER]))
        g = float(np.mean([out[f"{f}|{sname}|gain"]["cos_agg"] for f in CANCER]))
        dp = float(np.mean([out[f"{f}|{sname}|gain"]["pds"] - out[f"{f}|{sname}|none"]["pds"]
                            for f in CANCER]))
        summary[sname] = dict(cos_base=b, cos_gain=g, rel=100.0 * (g - b) / abs(b), d_pds=dp)
        print(f"{sname:<12}{b:>14.6f}{g:>15.6f}{100.0 * (g - b) / abs(b):>+9.3f}{dp:>+10.5f}")
    (WORK / "stack228.json").write_text(json.dumps(dict(rows=out, summary=summary), indent=2))
    print(f"-> {WORK / 'stack228.json'}")


if __name__ == "__main__":
    main()
