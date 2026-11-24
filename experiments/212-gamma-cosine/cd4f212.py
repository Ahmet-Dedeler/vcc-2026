"""212 — the same gamma sweep on `cd4f034`'s delivered profile (210's built arm).

`cd4f034` is `rc_callset` with `VCC_CD4_TARGET_F` 0.10 -> 0.34 and nothing else. The CD4
residual is added AFTER the operator and its fraction is re-solved from a ratio of median row
norms, so it interacts with gamma weakly but not at zero. 210 measured its delivered `A` at
**0.97970x** the incumbent's, which moves context A's rho cliff 1.03% further away
(`s_cliff` scales as `1/sqrt(0.9797)`), and the coordinator asked whether the gamma optimum
moves with it.

Identical to `plolo212.py` except `PROD["cd4_f"] = 0.34`. The `cd4_f = 0.10` column is the
known-answer row: it must reproduce `plolo212.json` bit for bit.

    uv run python cd4f212.py
"""
from __future__ import annotations

import gc
import json
import sys
import time
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
from plolo212 import DEV_KA, GRID, score_at, sha256_file  # noqa: E402
from plolo import FOLDS, PROD_OPS, Operator, build_fold, cd4_for  # noqa: E402
from score204 import DONORS, PROD, fold_sources, weighted_baseline  # noqa: E402

CD4_F = 0.34                # 210's `cd4f034`; the builder solves fraction 0.32207201067578145


def main() -> None:
    out_path = WORK / "cd4f212.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["gamma"]) for r in rows}
    PROD["cd4_f"] = CD4_F
    t0 = time.time()
    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        srcs = fold_sources(held, fold)
        cd4 = cd4_for(fold["targets"])
        op_path = Path(fold["op_path"] if FOLDS[held]["op_kind"] != "production"
                       else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        pre, _ = weighted_baseline(srcs, fold["targets"], "equal",
                                   {d: 1.0 for d in DONORS})
        dev = float(np.max(np.abs(pre - np.asarray(fold["pre"], dtype=np.float64))))
        assert dev == DEV_KA[held], f"`equal` drifted on {held}: {dev!r}"
        for gamma in GRID:
            if (held, gamma) in done:
                continue
            e = score_at(fold, pre, op, cd4, fold["scored"], truth, gamma)
            rows.append(dict(fold=held, gamma=gamma, cd4_f=CD4_F,
                             operator_sha256=sha256_file(op_path),
                             n_targets=len(fold["targets"]), **e))
            out_path.write_text(json.dumps(rows, indent=1))
            print(f"[{held}/g{gamma:.3f}] cos {e['a15_cos']:.7f} pds {e['a23_pds']:.6f} "
                  f"f {e['cd4_solved_fraction']:.6f} ({time.time()-t0:.0f}s)", flush=True)
        del op, srcs, cd4, pre
        gc.collect()

    folds = sorted(FOLDS)
    by = {(r["fold"], r["gamma"]): r for r in rows}
    base = {f: by[(f, 0.65)] for f in folds}

    def macro(g, key):
        return float(np.mean([by[(f, g)][key] for f in folds]))

    print(f"\n=== cd4_f {CD4_F} vs the incumbent's 0.10, macro over four folds ===")
    prim = {(r["fold"], r["gamma"]): r
            for r in json.loads((WORK / "plolo212.json").read_text())}

    def macro0(g, key):
        return float(np.mean([prim[(f, g)][key] for f in folds]))

    print("gamma |   cos(f=.34)  ratio |   cos(f=.10)  ratio | pds(.34)  d | pds(.10)  d")
    out = []
    for g in GRID:
        c34, c10 = macro(g, "a15_cos"), macro0(g, "a15_cos")
        p34, p10 = macro(g, "a23_pds"), macro0(g, "a23_pds")
        r34 = c34 / macro(0.65, "a15_cos")
        r10 = c10 / macro0(0.65, "a15_cos")
        out.append(dict(gamma=g, cos_f034=c34, ratio_f034=r34, cos_f010=c10,
                        ratio_f010=r10, pds_f034=p34, pds_f010=p10,
                        d_pds_f034=p34 - macro(0.65, "a23_pds"),
                        d_pds_f010=p10 - macro0(0.65, "a23_pds")))
        print(f"{g:5.3f} | {c34:.7f} {r34:6.4f} | {c10:.7f} {r10:6.4f} | "
              f"{p34:.6f} {p34-macro(0.65,'a23_pds'):+.6f} | "
              f"{p10:.6f} {p10-macro0(0.65,'a23_pds'):+.6f}")
    with open(WORK / "cd4f212_summary.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\n[cd4f] cd4_f 0.34 raises the absolute cosine by "
          f"{macro(0.65,'a15_cos')/macro0(0.65,'a15_cos')-1:+.2%} at gamma 0.65; the "
          f"gamma RATIO at 0.50 moves from {out[[o['gamma'] for o in out].index(0.50)]['ratio_f010']:.4f}"
          f" to {out[[o['gamma'] for o in out].index(0.50)]['ratio_f034']:.4f}")


if __name__ == "__main__":
    main()
