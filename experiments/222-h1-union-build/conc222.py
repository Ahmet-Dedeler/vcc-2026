"""222 INBOX-1 item 2 — is the union's gain carried by a handful of perturbations?

220 found this morning that both published ordering effects were carried by one to three
perturbations.  `PREREGISTRATION-ADDENDUM.md` section 5 pre-commits the reading:

    if the top three of the 77 carried targets carry more than half of `d pds` on the cancer
    folds, a positive board result is WEAKER evidence for 223's stack, because the board would
    be measuring a handful of perturbations rather than a source.

Both endpoints are per-target arrays out of 218's own `score_effects`, so nothing new is
defined here: `d pds` is the per-target PDS of the arm minus the baseline's, and `d cos` the
per-target cosine.  Leave-top-k-out is computed by DELETING those targets from the average,
which is the statistic the pre-commitment names.  The three cancer folds only -- `kolf2` is
excluded from every headline in 218 and here.

    uv run python conc222.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E218 = ROOT / "experiments/218-independent-donors"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(E218))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from plolo import (E062, FOLDS, OFFICIAL, PANEL, PROD_OPS, Operator,  # noqa: E402
                   build_fold, cd4_for, solve_fraction)
from resid218 import ALPHA, GAMMA, LAM, score_effects  # noqa: E402
from compose218 import merged  # noqa: E402

CD4_F, G = 0.34, 0.28
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")


def load_builder():
    path = ROOT / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("conc222_builder", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["conc222_builder"] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    MOD = load_builder()
    per_fold = {}
    for held in CANCER:
        fold = build_fold(held)
        targets = list(fold["targets"])
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(targets)
        union, _ = MOD.prepare_h1_union(targets, OFFICIAL, PANEL, max_rss_gib=40.0)
        carry = np.asarray(union.target_condition_count) > 0
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        whitened = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)
        frac_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
        frac_x, _ = MOD.solve_carrying_row_fraction(whitened, union, G)

        scored = {}
        for tag, residual, fx in (("base", None, 0.0), ("arm", union, frac_x)):
            final, _ = E062.add_post_whiten_residual(
                whitened, merged(cd4, frac_cd4, residual, fx), 1.0, fold["wcols"])
            final[:, fold["panel_positions"]] = 0.0
            scored[tag] = score_effects(final, fold, truth)
            del final
        d_pds = (np.asarray(scored["arm"]["_pds_per_target_a23"])
                 - np.asarray(scored["base"]["_pds_per_target_a23"]))
        d_cos = (np.asarray(scored["arm"]["_cos_per_target_a23"])
                 - np.asarray(scored["base"]["_cos_per_target_a23"]))
        per_fold[held] = dict(targets=targets, carry=carry.tolist(),
                              d_pds=d_pds.tolist(), d_cos=d_cos.tolist(),
                              base_cos=np.asarray(
                                  scored["base"]["_cos_per_target_a23"]).tolist(),
                              fold_d_pds=float(scored["arm"]["pds"]
                                               - scored["base"]["pds"]),
                              fold_rel_macro=float(scored["arm"]["a23_cos"]
                                                   / scored["base"]["a23_cos"] - 1))
        print(f"[{held}] {len(targets)} targets, {int(carry.sum())} carried; "
              f"fold d pds {per_fold[held]['fold_d_pds']:+.6f}, "
              f"rel macro cos {100*per_fold[held]['fold_rel_macro']:+.3f}%", flush=True)
        del op, whitened, union, cd4

    print("\n" + "=" * 104)
    print("CONCENTRATION — share of the mean d pds carried by the top-k targets, "
          "and the leave-top-k-out mean")
    print("=" * 104)
    rows = []
    for held in CANCER:
        f = per_fold[held]
        d_pds = np.asarray(f["d_pds"])
        d_cos = np.asarray(f["d_cos"])
        carry = np.asarray(f["carry"])
        order = np.argsort(-d_pds)
        total = float(d_pds.mean())
        print(f"\n{held}: mean d pds {total:+.6f} over {d_pds.size} targets "
              f"({int(carry.sum())} carried); "
              f"positive on {int((d_pds > 0).sum())}, negative on {int((d_pds < 0).sum())}")
        print(f"  top-5 by d pds: " + ", ".join(
            f"{f['targets'][i]}{'*' if carry[i] else ''} {d_pds[i]:+.4f}"
            for i in order[:5]))
        for k in (1, 3, 5, 10):
            top = order[:k]
            share = float(d_pds[top].sum() / d_pds.sum()) if d_pds.sum() != 0 else np.nan
            keep = np.ones(d_pds.size, dtype=bool)
            keep[top] = False
            print(f"  top-{k:<2d} share of total d pds {100*share:6.1f}%   "
                  f"leave-top-{k}-out mean d pds {float(d_pds[keep].mean()):+.6f}   "
                  f"(mean d cos {float(d_cos[keep].mean()):+.6f})")
            rows.append(dict(fold=held, k=k, share=share,
                             loo_mean_d_pds=float(d_pds[keep].mean()),
                             loo_mean_d_cos=float(d_cos[keep].mean()),
                             full_mean_d_pds=total))
        print(f"  carried rows carry {100*float(d_pds[carry].sum()/d_pds.sum()):5.1f}% of "
              f"the total d pds on {100*carry.mean():.1f}% of the targets")

    print("\n" + "=" * 104)
    print("POOLED over the three cancer folds")
    print("=" * 104)
    full = float(np.mean([per_fold[h]["fold_d_pds"] for h in CANCER]))
    print(f"  mean d pds, all targets                 {full:+.6f}")
    pooled = {}
    for k in (1, 3, 5, 10):
        vals = [r["loo_mean_d_pds"] for r in rows if r["k"] == k]
        shares = [r["share"] for r in rows if r["k"] == k]
        pooled[k] = dict(loo_mean_d_pds=float(np.mean(vals)),
                         mean_share=float(np.mean(shares)),
                         retained=float(np.mean(vals) / full) if full else np.nan)
        print(f"  leave-top-{k:<2d}-out mean d pds            "
              f"{np.mean(vals):+.6f}   ({100*np.mean(vals)/full:5.1f}% retained, "
              f"top-{k} mean share {100*np.mean(shares):5.1f}%)")
    verdict = pooled[3]["mean_share"] <= 0.5
    print(f"\nPRE-COMMITTED READING (addendum section 5): top-3 mean share "
          f"{100*pooled[3]['mean_share']:.1f}% "
          f"{'<= 50%  ->  the gain is DISTRIBUTED, the positive-board reading for 223 stands'
             if verdict else '> 50%  ->  CONCENTRATED, a positive board result is WEAKER evidence'}")

    (WORK / "conc222.json").write_text(json.dumps(
        dict(per_fold={h: {k: v for k, v in per_fold[h].items()} for h in CANCER},
             rows=rows, pooled={str(k): v for k, v in pooled.items()},
             full_mean_d_pds=full, top3_share_le_half=bool(verdict)), indent=2))
    print(f"\nwrote {WORK/'conc222.json'}")


if __name__ == "__main__":
    main()
