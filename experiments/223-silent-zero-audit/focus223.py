"""223 Phase 2b/3 — the two survivors, measured properly.

The broad sweep leaves exactly two sources positive on cosine AND ahead of both controls on
the three cancer folds: `dld1_lowmoi` (67/300) and `myllia_a549` (14/300).  This script:

  * sweeps a finer dose grid on both, on all four folds;
  * measures each ALONE on top of `cd4@0.34` as well as MARGINAL on top of
    `cd4@0.34 + h1pair@0.28`, so 222's contribution is never double-counted and the
    standalone value is available if 222 does not ship;
  * measures their UNION (their panel rows are almost disjoint);
  * repeats falsifier 1 with a SECOND independent derangement seed;
  * runs 218's / 207's leave-one-fold-out pricing: the dose is chosen on the other three
    folds and the price is read on the held one, so selection and reporting never overlap.

    uv run python focus223.py
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   solve_fraction)

sys.path.insert(0, str(HERE))
import sweep223  # noqa: E402  (registers dld1/altair/cmap in 218's source registry)
from sweep223 import audit, carrying_mask, row_split  # noqa: E402

sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
from resid218 import (ALPHA, GAMMA, LAM, score_effects, shuffled,  # noqa: E402
                      source_residual)
from compose218 import (merged, restrict_rows, solve_fraction_carrying,  # noqa: E402
                        union_rows)

CD4_F = 0.34
UNION_G = 0.28
UNION = ("gse295214_h1", "h1_2025")
SURV = tuple((__import__("os").environ.get("VCC223_SURV") or "dld1_lowmoi,myllia_a549").split(","))
DOSE = (0.05, 0.10, 0.17, 0.22, 0.28, 0.35, 0.45, 0.60)
SEED2 = 771103


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--out", default=str(WORK / "focus223.json"))
    args = ap.parse_args()

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in args.folds:
        fold = build_fold(held)
        targets = fold["targets"]
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(targets)
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        whitened = Operator(op_path).apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)

        def finish(res, fr):
            f = (whitened.copy() if res is None or fr == 0.0
                 else E062.add_post_whiten_residual(whitened, res, 1.0, fold["wcols"])[0])
            f[:, fold["panel_positions"]] = 0.0
            return f

        frac_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
        h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
        frac_u, _ = solve_fraction_carrying(whitened, h1, UNION_G)

        bases = {
            "cd4": merged(cd4, frac_cd4, None, 0.0),
            "base222": merged(cd4, frac_cd4, h1, frac_u),
        }
        base_e, base_cos_pt = {}, {}
        for bname, bres in bases.items():
            e = score_effects(finish(bres, 1.0), fold, truth)
            base_e[bname] = e
            base_cos_pt[bname] = np.asarray(e["_cos_per_target_a23"])
            if (held, f"BASE:{bname}") not in done:
                rows.append(dict(fold=held, arm=f"BASE:{bname}", kind="baseline",
                                 base=bname, **e))
            print(f"[{held}/BASE:{bname}] cos {e['cos']:.6f} pds {e['pds']:.6f}",
                  flush=True)

        srcs = {s: source_residual(s, targets, "gene_norm") for s in SURV if s != held}
        if len(srcs) == 2:
            srcs["dld1_x_myllia"] = union_rows(list(srcs.values()))
        variants = {}
        for sname, S in srcs.items():
            variants[sname] = S
            variants[f"{sname}|SHUF"] = shuffled(S)
            variants[f"{sname}|SHUF2"] = shuffled(S, seed=SEED2)
            variants[f"{sname}|CD4ROWS"] = restrict_rows(cd4, S)

        for bname, bres in bases.items():
            for vname, res in variants.items():
                cm = carrying_mask(res)
                for g in DOSE:
                    arm = f"{bname}+{vname}@{g:g}"
                    if (held, arm) in done:
                        continue
                    t0 = time.time()
                    fr, per_unit = solve_fraction_carrying(whitened, res, g)
                    assert fr > 0.0, f"vacuous arm {held}/{arm}"
                    m = merged(bres, 1.0, res, fr)
                    e = score_effects(finish(m, 1.0), fold, truth)
                    b = base_e[bname]
                    aud = audit(res, whitened, fr)
                    aud["per_unit"] = float(per_unit)
                    rows.append(dict(
                        fold=held, arm=arm, kind="arm", base=bname, source=vname, g=g,
                        n_targets=len(targets),
                        d_cos=e["cos"] - b["cos"],
                        rel_cos_pct=100.0 * (e["cos"] - b["cos"]) / max(abs(b["cos"]), 1e-300),
                        d_cos_agg=e["cos_agg"] - b["cos_agg"],
                        d_pds=e["pds"] - b["pds"], d_qw=e["a15_sgn_w"] - b["a15_sgn_w"],
                        rowsplit=row_split(e["_cos_per_target_a23"], base_cos_pt[bname], cm),
                        audit=aud, **e))
                    out_path.write_text(json.dumps(rows, indent=2))
                    print(f"[{held}/{arm}] frac {fr:.5f} carry {aud['carrying_rows']}"
                          f"/{aud['n_rows']} rel dcos {100.0 * (e['cos'] - b['cos']) / max(abs(b['cos']), 1e-300):+.3f}%"
                          f" dpds {e['pds'] - b['pds']:+.5f} ({time.time() - t0:.0f}s)",
                          flush=True)
        del whitened

    out_path.write_text(json.dumps(rows, indent=2))
    print(f"\n-> {out_path}")


if __name__ == "__main__":
    main()
