"""213 — is 212's 96.2% cosine-for-retrieval exchange LAW-INVARIANT?

212 6 closed `WHITEN_GAMMA` on a structural argument: along the gamma axis the operator converts
bulk cosine into retrieval at **96.2% of the board's own exchange rate** between `mse` and `pds`,

    d(overall mse*)/d(cos)  = +0.3655        (8.2's rate)
    d(overall pds )/d(cos)  = -0.3517        (what whitening charges, on the gamma axis)
                      sum   = +0.0138        (3.8%)

and it explicitly flagged that the ratio "need not be law-invariant" and should be re-checked
under any law that wins.  This measures it: the same gamma sweep, on the same four PLOLO folds,
under the incumbent law and under the family's cosine maximiser.

    uv run python gamma213.py
"""
from __future__ import annotations

import gc
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/204-repression-weighted-mixing"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402
plolo.WORK = WORK / "plolo"
import score204  # noqa: E402
score204.WORK = plolo.WORK.parent

from plolo import FOLDS, PROD_OPS, Operator, build_fold, cd4_for  # noqa: E402
import law213 as L  # noqa: E402
import plolo213 as S  # noqa: E402

GRID = [0.55, 0.60, 0.65, 0.70, 0.75]
LAWS = ["logc20", "logc40"]
FOLD_ORDER = ["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"]
COS0 = 0.077305107890            # chal213's own solved board level at the incumbent
RATE_PDS = 1.0 / 0.452 / 6.0 * 1.011
MSE_RATE_BOARD = 0.3655          # 8.2 / 212 3.4, the board's own rate


def main() -> None:
    laws = L.all_laws()
    out_path = WORK / "gamma213.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["law"], r["gamma"]) for r in rows}
    for held in FOLD_ORDER:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        cd4 = cd4_for(fold["targets"])
        op = Operator(Path(fold["op_path"] if FOLDS[held]["op_kind"] != "production"
                           else PROD_OPS / "whitener_union_pruned.npz"))
        for key in LAWS:
            srcs = S.law_sources(held, fold, laws[key])
            pre, _ = score204.weighted_baseline(srcs, fold["targets"], "equal",
                                                {n: 1.0 for n, _ in srcs})
            for g in GRID:
                if (held, key, g) in done:
                    continue
                e = S.score_law(fold, pre, op, cd4, fold["scored"], truth, laws[key],
                                gamma=g)
                rows.append(dict(fold=held, law=key, gamma=g,
                                 a15_cos=e["a15_cos"], a15_pds=e["a15_pds"],
                                 a23_cos=e["a23_cos"], a23_pds=e["a23_pds"]))
                out_path.write_text(json.dumps(rows, indent=1))
                print(f"[{held}/{key}/g{g}] cos {e['a15_cos']:.7f} pds {e['a15_pds']:.6f}",
                      flush=True)
            del srcs, pre
            gc.collect()
        del op, cd4
        gc.collect()

    def macro(law, g, key):
        v = {r["fold"]: r[key] for r in rows if r["law"] == law and r["gamma"] == g}
        return float(np.mean([v[f] for f in FOLD_ORDER]))

    print("\n=== 212 6's exchange, re-measured under each law (gamma 0.60 -> 0.70) ===")
    res = {}
    for law in LAWS:
        c_lo, c_hi = macro(law, 0.70, "a15_cos"), macro(law, 0.60, "a15_cos")
        p_lo, p_hi = macro(law, 0.70, "a15_pds"), macro(law, 0.60, "a15_pds")
        # board cosine scale: the law's own PLOLO cosine at gamma 0.65 maps to COS0 * ratio
        c65 = macro(law, 0.65, "a15_cos")
        c65_inc = macro("logc20", 0.65, "a15_cos")
        board = COS0 * c65 / c65_inc
        dcos_board = board * (c_hi - c_lo) / c65
        dpds_overall = (p_hi - p_lo) * RATE_PDS
        rate = dpds_overall / dcos_board
        res[law] = dict(board_cos_at_065=board, d_board_cos=dcos_board,
                        d_overall_pds=dpds_overall, d_overall_pds_per_cos=rate,
                        board_mse_rate=MSE_RATE_BOARD,
                        sum_rate=MSE_RATE_BOARD + rate,
                        cancellation_pct=100.0 * (1.0 + rate / MSE_RATE_BOARD))
        print(f"  {law:<8} board cos {board:.6f} | d(overall pds)/d(cos) = {rate:+.4f} "
              f"vs the board's mse rate +{MSE_RATE_BOARD} -> sum {MSE_RATE_BOARD + rate:+.4f} "
              f"({res[law]['cancellation_pct']:.1f}% of the mse leg survives)")
    (WORK / "gamma213_rate.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
