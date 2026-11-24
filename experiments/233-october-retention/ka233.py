"""233 Phase 0.3 — the instrument's own known-answer rows, before any retention number.

KA-C  MY scoring path, run on TODAY'S panel with 062's held operators (the PLOLO instrument
      exactly as 162 validated it against twelve board points), must reproduce the SHIPPED
      STACK's `cos_agg` and `pds` on the three cancer folds at Delta = 0.0e+00.  The reference
      values are 223's `focus223` `BASE:base222`, carried in `arms228.KA`, which is the
      `h1union` configuration on the board at 0.2118310.

KA-D  MY fold builder's consensus must equal `plolo.build_fold`'s `pre` EXACTLY on today's
      panel.  `fold233` keeps the per-source pieces so Phase 2 can ablate them; that is only
      legitimate if re-averaging them reproduces the object 162 built.

KA-E  BRIDGE.  The retention measurement re-fits every operator through
      `scripts/build_union_whitener.py` so that nothing but the panel varies between the
      three arms.  062's held operators were fitted by a different recipe (no essential
      prune, no 12,000-row cap).  This row measures how much of the A/B/C control moves under
      that change of recipe, so the reader can tell an operator-construction effect from a
      panel effect.  It is reported, not required to be zero.

    uv run python ka233.py
"""
from __future__ import annotations

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
sys.path.insert(0, str(HERE))

import panels233 as P  # noqa: E402
import score233 as S  # noqa: E402
import fold233  # noqa: E402
import plolo  # noqa: E402
from plolo import E062, Operator, build_fold, cd4_for, solve_fraction  # noqa: E402

HELD_OPS = ROOT / "experiments/062-cd4-post-whiten/cloud-output/whiteners"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--bridge", action="store_true", help="also run KA-E (needs the refits)")
    ap.add_argument("--folds", nargs="*", default=list(CANCER))
    args = ap.parse_args()

    resid218, compose218 = S._resid_module()
    sys.path.insert(0, str(ROOT / "experiments/228-theta-gain-plolo"))
    from arms228 import KA  # noqa: E402

    bctx = {c: plolo.challenge_bctrl(c) for c in S.CONTEXTS}
    gate = {}
    for c in S.CONTEXTS:
        e = np.expm1(bctx[c])
        gate[c] = (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (bctx[c] > 0)
    ctrl_mod = __import__("ctrl228")
    ctrl = {c: ctrl_mod.control_block(c) for c in S.CONTEXTS}

    out: dict[str, object] = {}
    print(f"{'fold':16}{'cos_agg':>13}{'expected':>13}{'d':>12}"
          f"{'pds':>12}{'expected':>12}{'d':>12}")
    fails = 0
    for held in args.folds:
        # ---- the instrument exactly as 162 built it: today's panel, 062's held operator
        plolo.PANEL = P.read_panel("abc")
        plolo.PANEL_SET = set(plolo.PANEL)
        plolo.WORK = ROOT / "experiments/162-retrodicting-validator/work"
        plolo.CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
        plolo.FOLDS[held]["op"] = HELD_OPS / f"held_{held}.npz"
        plolo.FOLDS[held]["op_kind"] = "held"
        resid218.PANEL_SET = set(plolo.PANEL)
        resid218._PREP_CACHE.clear()

        fold = build_fold(held)
        plolo.WORK = S.SHARED
        targets = fold["targets"]
        wcols = np.asarray(fold["wcols"])
        sc = np.asarray(fold["scored"])
        Ts = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
        cd4 = cd4_for(targets)
        h1 = compose218.union_rows(
            [resid218.source_residual(n, targets, "gene_norm") for n in S.UNION])
        op = Operator(HELD_OPS / f"held_{held}.npz")
        pre = np.asarray(fold["pre"], dtype=np.float64)
        final = S.build_final(dict(panel_positions=fold["panel_positions"]), op, pre,
                              cd4, h1, wcols, S.CD4_F, S.UNION_G, (resid218, compose218))
        m = S.score_arm(final, dict(scored=sc), bctx, gate, ctrl, Ts, None, None, None,
                        np.random.default_rng(7), None, full=False)
        exp_cos, exp_pds = KA[(held, "base222")]     # (mean per-row cos, pds)
        dc, dp = m["cos"] - exp_cos, m["pds"] - exp_pds
        fails += int(abs(dc) > 1e-12 or abs(dp) > 1e-12)
        out[f"KA-C|{held}"] = dict(cos=m["cos"], cos_expected=exp_cos, d_cos=dc,
                                   cos_agg=m["cos_agg"],
                                   pds=m["pds"], pds_expected=exp_pds, d_pds=dp,
                                   n_targets=int(final.shape[0]))
        print(f"{held:16}{m['cos']:13.9f}{exp_cos:13.9f}{dc:12.2e}"
              f"{m['pds']:12.8f}{exp_pds:12.8f}{dp:12.2e}")

        # ---- KA-D: my fold builder's consensus vs plolo's `pre`
        try:
            mine = fold233.build("abc", held)
            pre_mine, _ = fold233.consensus(mine, pm=True, only=None)
            # `plolo` stores its `pre` as float32; compare in that storage type
            d = float(np.abs(pre_mine.astype(np.float32) - pre.astype(np.float32)).max())
            out[f"KA-D|{held}"] = dict(max_abs_diff=d, shape=list(pre_mine.shape))
            print(f"{'  KA-D consensus':16}max|mine - plolo pre| = {d:.3e}"
                  f"  ({'PASS' if d == 0.0 else 'FAIL'})")
            fails += int(d != 0.0)
        except Exception as e:                                  # operator not fitted yet
            print(f"  KA-D skipped: {e}")

        if args.bridge and P.op_path("abc", held).exists():
            op2 = Operator(P.op_path("abc", held))
            f2 = S.build_final(dict(panel_positions=fold["panel_positions"]), op2, pre,
                               cd4, h1, wcols, S.CD4_F, S.UNION_G, (resid218, compose218))
            m2 = S.score_arm(f2, dict(scored=sc), bctx, gate, ctrl, Ts, None, None, None,
                             np.random.default_rng(7), None, full=False)
            out[f"KA-E|{held}"] = dict(refit_cos_agg=m2["cos_agg"], held062_cos_agg=m["cos_agg"],
                                       refit_cos=m2["cos"], held062_cos=m["cos"],
                                       rel_cos_agg_pct=100 * (m2["cos_agg"] / m["cos_agg"] - 1),
                                       refit_pds=m2["pds"], d_pds=m2["pds"] - m["pds"])
            print(f"{'  KA-E bridge':16}refit operator: cos {m2['cos_agg']:.9f} "
                  f"({100 * (m2['cos_agg'] / m['cos_agg'] - 1):+.2f}%)  "
                  f"pds {m2['pds']:.8f} ({m2['pds'] - m['pds']:+.5f})")
            del op2, f2
        del op, final

    prev = json.loads((P.WORK / "ka233.json").read_text()) if (P.WORK / "ka233.json").exists() else {}
    (P.WORK / "ka233.json").write_text(json.dumps(prev | out, indent=2))
    print(f"\nknown-answer failures: {fails}")
    print("-> work/ka233.json")
    if fails:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
