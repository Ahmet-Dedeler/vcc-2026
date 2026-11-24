"""243 Phase 1 — the PLOLO side of the three board direction steps, every way it can be
defined, read from the artifacts that already hold it.

  step        before -> after                    artifact (per fold, both tags)
  a14->a15    a14 (fraction 0.05)  -> f 0.10     experiments/207-cd4-dose/work/sweep207.json
  cd4f034     f 0.10 -> f 0.34                   experiments/216-cd4-hard/work/sweep216.json
  h1union     f 0.34 -> f 0.34 + union g 0.28    experiments/218-independent-donors/work/
                                                 compose218.json (a15+a23 tags) and 222's
                                                 repro222.json (a23 tag, production union,
                                                 reproduces 218 to 0.0e+00)

Axes: functional {cos_agg, macro} x tag {a23, a15} x fold set {4, 3 cancer} x aggregation
{pooled = ratio of fold-means (216), perfold = mean of per-fold ratios (222)}.

Also `d pds` on both fold sets and both tags.  Board-side numbers are attached in rule243.py.

    uv run python matrix243.py
"""
from __future__ import annotations

import json

import numpy as np

from lib243 import (CANCER, FOLDS, FUNCTIONALS, PDS_KEYS, WORK, compose218, d_pds, rel_gain,
                    repro222, sweep207, sweep216)


def main() -> None:
    s207, s216, c218 = sweep207(), sweep216(), compose218()
    r222 = repro222()
    a222 = {(r["fold"], r["family"], r["g"]): r for r in r222["arms"]}

    # --- known-answer consistency across the artifacts (same object, three writers) --------
    checks = []
    for f in FOLDS:
        # 207's f0.1 row == 216's f0.1 row (216 §4.1: "every one of 207's rows reproduces")
        checks.append(("207 f0.1 == 216 f0.1 cos_agg", f,
                       abs(s207[(f, "f0.1")]["cos_agg"] - s216[(f, "f0.1")]["cos_agg"])))
        checks.append(("207 f0.34 == 216 f0.34 cos_agg", f,
                       abs(s207[(f, "f0.34")]["cos_agg"] - s216[(f, "f0.34")]["cos_agg"])))
        # 218's base row == 216's f0.34 row
        checks.append(("218 cd40.34+none0 == 216 f0.34 cos_agg", f,
                       abs(c218[(f, "cd40.34+none0")]["cos_agg"] - s216[(f, "f0.34")]["cos_agg"])))
        # 222's production union == 218's union (222 KA-3b)
        u = a222[(f, "h1PAIRCARRY", 0.28)]
        checks.append(("222 union == 218 union cos_agg", f,
                       abs(u["cos_agg"] - c218[(f, "cd40.34+h1PAIRCARRY0.28")]["cos_agg"])))
        checks.append(("222 union == 218 union a23_cos", f,
                       abs(u["a23_cos"] - c218[(f, "cd40.34+h1PAIRCARRY0.28")]["a23_cos"])))
        checks.append(("222 base == 218 base cos_agg", f,
                       abs(u["base_cos_agg"] - c218[(f, "cd40.34+none0")]["cos_agg"])))
    worst = max(d for _, _, d in checks)
    print(f"artifact cross-checks: {len(checks)} rows, worst |d| = {worst:.2e}  "
          f"{'PASS' if worst < 1e-9 else 'FAIL'}")
    if worst >= 1e-9:
        for n, f, d in checks:
            if d >= 1e-9:
                print("   ", n, f, d)
        raise SystemExit(1)

    steps = {
        "a14->a15": ({f: s207[(f, "a14")] for f in FOLDS}, {f: s207[(f, "f0.1")] for f in FOLDS}),
        "cd4f034": ({f: s216[(f, "f0.1")] for f in FOLDS}, {f: s216[(f, "f0.34")] for f in FOLDS}),
        "h1union": ({f: c218[(f, "cd40.34+none0")] for f in FOLDS},
                    {f: c218[(f, "cd40.34+h1PAIRCARRY0.28")] for f in FOLDS}),
    }
    out = {"cells": [], "per_fold": {}}
    print(f"\n{'step':10s} {'functional':10s} {'tag':4s} {'folds':6s} {'agg':8s} {'rel gain %':>11s}")
    for step, (before, after) in steps.items():
        out["per_fold"][step] = {}
        for (func, tag), key in FUNCTIONALS.items():
            out["per_fold"][step][f"{func}|{tag}"] = {
                f: dict(before=before[f][key], after=after[f][key],
                        rel=after[f][key] / before[f][key] - 1.0) for f in FOLDS}
            for fs_name, fs in (("4", FOLDS), ("3c", CANCER)):
                for how in ("pooled", "perfold"):
                    rel = rel_gain(before, after, key, fs, how)
                    out["cells"].append(dict(step=step, functional=func, tag=tag, folds=fs_name,
                                             agg=how, rel=rel))
                    print(f"{step:10s} {func:10s} {tag:4s} {fs_name:6s} {how:8s} {100*rel:>+11.4f}")
        for tag, key in PDS_KEYS.items():
            for fs_name, fs in (("4", FOLDS), ("3c", CANCER)):
                dp = d_pds(before, after, key, fs)
                out["cells"].append(dict(step=step, functional="pds", tag=tag, folds=fs_name,
                                         agg="mean", rel=dp))
                print(f"{step:10s} {'d pds':10s} {tag:4s} {fs_name:6s} {'mean':8s} {dp:>+11.6f}")
        print()

    # per-fold table for the RESULTS
    for step in steps:
        print(f"per fold, {step}, rel cos_agg (a23) / rel macro (a23) / d pds (a23):")
        for f in FOLDS:
            a = out["per_fold"][step]["cos_agg|a23"][f]["rel"]
            m = out["per_fold"][step]["macro|a23"][f]["rel"]
            dp = steps[step][1][f]["pds"] - steps[step][0][f]["pds"]
            print(f"   {f:14s} {100*a:+8.3f}%  {100*m:+8.3f}%  {dp:+.6f}")
    (WORK / "matrix243.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {WORK/'matrix243.json'}")


if __name__ == "__main__":
    main()
