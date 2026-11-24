"""252 known-answer row 1a -- 218's `compose218` headline rows through the guarded harness,
called EXACTLY as `compose218.py` calls it: `build_fold` / `cd4_for` / `Operator` /
`solve_fraction` for CD4 at 0.34 (and 0.10), `resid218.source_residual` + `compose218.union_rows`
for the H1 pair, `compose218.solve_fraction_carrying` for its dose, `compose218.merged`, ONE
`E062.add_post_whiten_residual`, `resid218.score_effects`.  Four folds; the arms

    cd40.34+none0   cd40.34+h1PAIRCARRY0.28   cd40.34+h1PAIRCARRY_SHUF0.28   cd40.1+none0

against `work/compose218.json` (per fold: cos_agg, a23_cos, pds, a15_sgn_w, frac_cd4, frac_x,
preserve_scale) and the headline on the three cancer folds against `work/diag218.json`:
rel_macro_nokolf +0.07680281456197664, d_pds_nokolf +0.004886737506532087.

    /usr/bin/time -l uv run python ka218_252.py        -> work/ka218.json ; exit 1 on any |d| > 0
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
E218 = ROOT / "experiments/218-independent-donors"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(E218))

import plolo  # noqa: E402
import compose218 as C  # noqa: E402  (points plolo.WORK at 218's caches)
import resid218 as R  # noqa: E402

plolo.WORK = WORK / "plolo"            # 252's byte copies (sha-identical to 218's)
from plolo import E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for, solve_fraction  # noqa: E402

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
KEYS = ("cos_agg", "a23_cos", "pds", "a15_sgn_w", "a15_cos")
HEAD_MACRO = 0.07680281456197664
HEAD_PDS = 0.004886737506532087


def main() -> None:
    t0 = time.time()
    ref = {(r["fold"], r["arm"]): r
           for r in json.loads((E218 / "work/compose218.json").read_text())}
    ref_macro = json.loads((E218 / "work/diag218.json").read_text())["macro"]
    rows, worst = [], 0.0
    for held in FOLDS:
        fold = build_fold(held)
        targets = fold["targets"]
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(targets)
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        whitened = op.apply(pre * R.ALPHA, fold["wcols"], R.GAMMA, R.LAM)
        op_rec = dict(name=op.name, sha256=op.sha256, axis_checked=op.axis_checked)
        del op
        gc.collect()
        pair = C.union_rows([R.source_residual(n, targets, "gene_norm") for n in C.SECOND_CARRY])
        seconds = {"h1PAIRCARRY": pair, "h1PAIRCARRY_SHUF": R.shuffled(pair)}
        for cf, arms in ((0.34, (("none", 0.0), ("h1PAIRCARRY", 0.28), ("h1PAIRCARRY_SHUF", 0.28))),
                         (0.10, (("none", 0.0),))):
            frac_cd4, _ = solve_fraction(whitened, cd4, cf)
            for sname, g in arms:
                arm = f"cd4{cf:g}+{sname}{g:g}"
                sres = seconds.get(sname)
                if sres is None or g == 0.0:
                    frac_x = 0.0
                else:
                    frac_x, _ = C.solve_fraction_carrying(whitened, sres, g)
                    assert frac_x > 0.0, f"VACUOUS second source {sname} on {held}"
                res = C.merged(cd4, frac_cd4, sres, frac_x)
                final, preserve = E062.add_post_whiten_residual(whitened, res, 1.0, fold["wcols"])
                final[:, fold["panel_positions"]] = 0.0
                e = R.score_effects(final, fold, truth)
                r = ref[(held, arm)]
                d = {k: abs(float(e[k]) - float(r[k])) for k in KEYS}
                d["frac_cd4"] = abs(float(frac_cd4) - float(r["audit"]["frac_cd4"]))
                d["frac_x"] = abs(float(frac_x) - float(r["audit"]["frac_x"]))
                d["preserve_scale"] = abs(float(preserve["post_add_global_preserve_scale"])
                                          - float(r["audit"]["preserve_scale"]))
                m = max(d.values())
                worst = max(worst, m)
                rows.append(dict(fold=held, arm=arm, cd4_f=cf, g=g, frac_cd4=float(frac_cd4),
                                 frac_x=float(frac_x), cos_agg=float(e["cos_agg"]),
                                 a23_cos=float(e["a23_cos"]), pds=float(e["pds"]),
                                 ref_cos_agg=r["cos_agg"], ref_a23_cos=r["a23_cos"], ref_pds=r["pds"],
                                 max_abs_diff=m, diffs=d, operator=op_rec,
                                 cache=fold["audit"].get("cache")))
                print(f"[{held}/{arm}] fx {frac_x:.6f} cos_agg {e['cos_agg']:.15g} "
                      f"a23_cos {e['a23_cos']:.15g} pds {e['pds']:.15g}  vs 218: max|d| {m:.1e}",
                      flush=True)
                del final
        del whitened, pre, cd4, fold, truth, seconds, pair
        gc.collect()

    # the headline: three cancer folds, macro cosine and pds relative to cd40.34+none0
    by = {(r["fold"], r["arm"]): r for r in rows}
    head = {}
    for name in ("h1PAIRCARRY", "h1PAIRCARRY_SHUF"):
        rm, dp = [], []
        for h in CANCER:
            b, a = by[(h, "cd40.34+none0")], by[(h, f"cd40.34+{name}0.28")]
            rm.append(a["a23_cos"] / b["a23_cos"] - 1.0)
            dp.append(a["pds"] - b["pds"])
        head[name] = dict(rel_macro_nokolf=float(np.mean(rm)), d_pds_nokolf=float(np.mean(dp)),
                          ref_rel_macro_nokolf=ref_macro[f"{name}@0.28"]["rel_macro_nokolf"],
                          ref_d_pds_nokolf=ref_macro[f"{name}@0.28"]["d_pds_nokolf"])
        head[name]["d_macro"] = abs(head[name]["rel_macro_nokolf"] - head[name]["ref_rel_macro_nokolf"])
        head[name]["d_pds"] = abs(head[name]["d_pds_nokolf"] - head[name]["ref_d_pds_nokolf"])
        worst = max(worst, head[name]["d_macro"], head[name]["d_pds"])
    hm = head["h1PAIRCARRY"]
    mandate_ok = (abs(hm["rel_macro_nokolf"] - HEAD_MACRO) == 0.0
                  and abs(hm["d_pds_nokolf"] - HEAD_PDS) == 0.0)
    print(f"HEADLINE cd40.34+h1PAIRCARRY0.28, 3 cancer folds: macro {100 * hm['rel_macro_nokolf']:+.4f}% "
          f"(mandate +7.6803%) pds {hm['d_pds_nokolf']:+.6f} (mandate +0.004887); "
          f"|d| vs diag218 {hm['d_macro']:.1e} / {hm['d_pds']:.1e}; "
          f"shuffled control macro {100 * head['h1PAIRCARRY_SHUF']['rel_macro_nokolf']:+.4f}%")
    out = dict(rows=rows, headline=head, headline_reproduces_mandate_digits=bool(mandate_ok),
               worst_abs_diff=worst, n_rows=len(rows), seconds=round(time.time() - t0, 1),
               peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3)
    (WORK / "ka218.json").write_text(json.dumps(out, indent=1))
    print(f"KA-218: {len(rows)} rows, worst |d| {worst:.1e}, {out['seconds']}s, "
          f"peak {out['peak_rss_gib']:.2f} GiB -> {'PASS' if worst == 0.0 else 'FAIL'}")
    if worst > 0.0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
