"""238 — WHICH TARGETS carry the residual gain on the essential-free panels, by name.

`six238.py`'s leave-top-3-out reports HOW MUCH of the gain three rows carry; this reports WHICH
three, because on a panel where the whole asset is 9 carrying rows the identity of those rows is
the finding.  Same construction as `six236.leave_top_out` (per-row leave-one-out influence on the
GAIN, `pds` dropping the row from the decoy set too), with the target names kept.

    uv run python influence238.py --panels efp efr pess
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
sys.path.insert(0, str(HERE))

import bind238  # noqa: E402,F401
import panels238 as P  # noqa: E402
import score233 as S  # noqa: E402
import fold233  # noqa: E402
import plolo  # noqa: E402
from plolo import E062, Operator, deliver  # noqa: E402
from six236 import cos_from, pds_from, prep_rows  # noqa: E402

ESSENTIAL = {ln.strip() for ln in
             (P.ROOT / "data/donors/k562_essential/perts.txt").read_text().splitlines()
             if ln.strip()}
SOURCES = ("nadig_hepg2", "nadig_jurkat", "rpe1")
FOLDS = ("orion_hct116", "orion_hek293t")
DOSE = 0.34


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="*", default=["efp", "efr", "pess"])
    args = ap.parse_args()
    resid218, compose218 = S._resid_module()
    bctx = {c: plolo.challenge_bctrl(c) for c in S.CONTEXTS}
    gate = {}
    for c in S.CONTEXTS:
        e = np.expm1(bctx[c])
        gate[c] = (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (bctx[c] > 0)

    out: dict = {}
    for panel in args.panels:
        S.bind_panel(panel, resid218)
        for held in FOLDS:
            fold = fold233.build(panel, held)
            sc = np.asarray(fold["scored"])
            wcols = np.asarray(fold["wcols"])
            Ts = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
            targets = fold["targets"]
            cd4 = S.cd4_prepared(panel, targets)
            h1 = compose218.union_rows(
                [resid218.source_residual(n, targets, "gene_norm") for n in S.UNION])
            op = Operator(fold["op_path"])
            pre = fold233.consensus(fold, pm=True, only=None)[0]
            wh = op.apply(pre * S.ALPHA, wcols, S.GAMMA, S.LAM)
            f_cd4, _ = plolo.solve_fraction(wh, cd4, S.CD4_F)
            f_u, _ = compose218.solve_fraction_carrying(wh, h1, S.UNION_G)
            base_resid = compose218.merged(cd4, f_cd4, h1, f_u)
            res = compose218.union_rows(
                [resid218.source_residual(n, targets, "gene_norm") for n in SOURCES])
            fx, _ = compose218.solve_fraction_carrying(wh, res, DOSE)

            def final_of(add):
                f, _ = E062.add_post_whiten_residual(wh, add, 1.0, wcols)
                f[:, fold["panel_positions"]] = 0.0
                return f

            fb = final_of(base_resid)
            fa = final_of(compose218.merged(base_resid, 1.0, res, fx))
            n = fb.shape[0]
            idx = np.arange(n)
            inf_cos = np.zeros(n)
            inf_pds = np.zeros(n)
            gain_cos = gain_pds = 0.0
            for c in S.CONTEXTS:
                prb = prep_rows(deliver(fb, bctx[c], gate[c], "minnorm")[:, sc], Ts)
                pra = prep_rows(deliver(fa, bctx[c], gate[c], "minnorm")[:, sc], Ts)
                gc = cos_from(pra, idx) - cos_from(prb, idx)
                gp = pds_from(pra, idx) - pds_from(prb, idx)
                gain_cos += gc / 3.0
                gain_pds += gp / 3.0
                for i in range(n):
                    m = np.delete(idx, i)
                    inf_cos[i] += (gc - (cos_from(pra, m) - cos_from(prb, m))) / 3.0
                    inf_pds[i] += (gp - (pds_from(pra, m) - pds_from(prb, m))) / 3.0
                del prb, pra
            carry = np.asarray(res.target_condition_count) > 0
            order = np.argsort(-inf_cos)
            # `share` is a row's leave-one-out influence as a fraction of the ARM'S ACTUAL
            # COSINE GAIN.  Influences do not sum to the gain (the aggregate cosine's
            # denominator moves when a row is dropped), so they are not a partition -- shares
            # can exceed 100% and need not add to it.  The honest summary of concentration is
            # the leave-top-3-out retention `six238.py` reports; this names the rows.
            tot = gain_cos
            rows = [dict(target=targets[i], carrying=bool(carry[i]),
                         essential=targets[i] in ESSENTIAL,
                         inf_cos=float(inf_cos[i]), inf_pds=float(inf_pds[i]),
                         share=float(inf_cos[i] / tot) if tot else float("nan"))
                    for i in order[:8]]
            out[f"{panel}|{held}"] = dict(
                n_targets=n, n_carrying=int(carry.sum()), gain_cos=float(gain_cos),
                gain_pds=float(gain_pds),
                top1=str(targets[order[0]]),
                top1_share=float(inf_cos[order[0]] / tot) if tot else float("nan"),
                top=rows)
            print(f"\n[{panel}/{held}]  {int(carry.sum())} carrying rows of {n}; "
                  f"cosine gain {gain_cos:+.6f}, pds gain {gain_pds:+.6f}; "
                  f"top row {targets[order[0]]} = "
                  f"{100 * inf_cos[order[0]] / tot:.0f}% of the gain")
            print(f"   {'target':14}{'carrying':>9}{'essential':>11}{'inf_cos':>12}"
                  f"{'share':>9}{'inf_pds':>11}")
            for r in rows:
                print(f"   {r['target']:14}{str(r['carrying']):>9}{str(r['essential']):>11}"
                      f"{r['inf_cos']:12.6f}{100 * r['share']:8.1f}%{r['inf_pds']:11.6f}")
            del op, wh, fb, fa
    (P.WORK / "influence238.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> {P.WORK / 'influence238.json'}")


if __name__ == "__main__":
    main()
