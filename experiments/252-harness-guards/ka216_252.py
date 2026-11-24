"""252 known-answer row 4 -- 216's `f = 0.34` and `f = 0.10` rows through the guarded harness,
called EXACTLY as `sweep216.py` calls it (its own `score_effects`, its own ALPHA/GAMMA/LAM,
`plolo.build_fold` / `Operator` / `solve_fraction` / `E062.add_post_whiten_residual`), on all
four folds, against `experiments/216-cd4-hard/work/sweep216.json`.

    /usr/bin/time -l uv run python ka216_252.py        -> work/ka216.json ; exit 1 on any |d| > 0
"""
from __future__ import annotations

import gc
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E216 = ROOT / "experiments/216-cd4-hard"
sys.path.insert(0, str(E216))

import sweep216 as S  # noqa: E402  (imports plolo from 162 and points plolo.WORK at 216's caches)
import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"            # 252's byte copies (sha-identical to 216's)
from plolo import E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for, solve_fraction  # noqa: E402

ARMS = (("f0.34", 0.34), ("f0.1", 0.10))
KEYS = ("cos_agg", "pds", "cos", "a15_cos", "a15_pds", "sgn_w", "sse")


def main() -> None:
    t0 = time.time()
    ref = {(r["fold"], r["arm"]): r
           for r in json.loads((E216 / "work/sweep216.json").read_text())}
    rows, worst = [], 0.0
    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        whitened = op.apply(pre * S.ALPHA, fold["wcols"], S.GAMMA, S.LAM)
        op_rec = dict(name=op.name, sha256=op.sha256, axis_checked=op.axis_checked)
        del op
        gc.collect()
        for name, f_val in ARMS:
            frac, per_unit = solve_fraction(whitened, cd4, f_val)
            final, preserve = E062.add_post_whiten_residual(whitened, cd4, frac, fold["wcols"])
            final[:, fold["panel_positions"]] = 0.0
            e = S.score_effects(final, fold, truth)
            r = ref[(held, name)]
            d = {k: abs(float(e[k]) - float(r[k])) for k in KEYS}
            d["fraction"] = abs(float(frac) - float(r["fraction"]))
            d["f_per_unit"] = abs(float(per_unit) - float(r["f_per_unit"]))
            d["preserve_scale"] = abs(float(preserve["post_add_global_preserve_scale"])
                                      - float(r["preserve_scale"]))
            m = max(d.values())
            worst = max(worst, m)
            rows.append(dict(fold=held, arm=name, f=f_val, fraction=float(frac),
                             cos_agg=float(e["cos_agg"]), pds=float(e["pds"]),
                             a15_cos=float(e["a15_cos"]),
                             ref_cos_agg=r["cos_agg"], ref_pds=r["pds"], ref_fraction=r["fraction"],
                             max_abs_diff=m, diffs=d, operator=op_rec,
                             cache=fold["audit"].get("cache")))
            print(f"[{held}/{name}] frac {frac:.15g} cos_agg {e['cos_agg']:.15g} "
                  f"pds {e['pds']:.15g}  vs 216: max|d| {m:.1e}", flush=True)
            del final
        del whitened, pre, cd4, fold, truth
        gc.collect()
    out = dict(rows=rows, worst_abs_diff=worst, n_rows=len(rows),
               seconds=round(time.time() - t0, 1),
               peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3)
    (WORK / "ka216.json").write_text(json.dumps(out, indent=1))
    print(f"KA-216: {len(rows)} rows, worst |d| {worst:.1e}, {out['seconds']}s, "
          f"peak {out['peak_rss_gib']:.2f} GiB -> {'PASS' if worst == 0.0 else 'FAIL'}")
    if worst > 0.0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
