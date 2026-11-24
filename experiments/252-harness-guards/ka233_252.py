"""252 known-answer row 3 -- 233's BASE rows on abc / real / pess through the guarded harness,
by running `score233.main()` ITSELF (its own argv: `--panels abc real pess --arms BASE`), with
one thing named up front:

    plolo.SIZING_DEFAULT = "legacy_all_rows"

because `score233.build_final` reads production's ALL-ROWS per-unit of the H1 union PURELY AS
A DIAGNOSTIC (`sizing.h1_per_unit_all_rows`, recorded 0.0 on all 81 of 233's rows: the union
covers 77/300 on every panel) and under the strict default that read raises `SilentZeroError`.
The APPLIED fractions never touch the zero: CD4 is sized at > 50% coverage on all three panels
(per_unit > 0, arithmetic bit-identical) and the union by 218's carrying rule.  Every
`[plolo] LEGACY SIZING` line in the log is one of those diagnostic reads.

Reads 233's fold caches / re-fit operators / CD4 banks; writes ONLY `work/score233-BASE.json`
and `work/ka233.json`.

    /usr/bin/time -l uv run python ka233_252.py        ; exit 1 on any |d| > 0
"""
from __future__ import annotations

import json
import math
import resource
import runpy
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E233 = ROOT / "experiments/233-october-retention"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.SIZING_DEFAULT = "legacy_all_rows"          # named, once, for the diagnostic reads only


def leaves(v, prefix=""):
    if isinstance(v, dict):
        out = {}
        for k, x in v.items():
            out.update(leaves(x, f"{prefix}{k}."))
        return out
    if isinstance(v, list):
        out = {}
        for i, x in enumerate(v):
            out.update(leaves(x, f"{prefix}{i}."))
        return out
    return {prefix.rstrip("."): v}


def main() -> None:
    t0 = time.time()
    out_path = WORK / "score233-BASE.json"
    if out_path.exists():
        out_path.unlink()
    sys.argv = ["score233.py", "--panels", "abc", "real", "pess", "--arms", "BASE",
                "--out", str(out_path)]
    sys.path.insert(0, str(E233))
    runpy.run_path(str(E233 / "score233.py"), run_name="__main__")

    ref = {(r["panel"], r["fold"]): r
           for r in json.loads((E233 / "work/score233.json").read_text()) if r["arm"] == "BASE"}
    got = {(r["panel"], r["fold"]): r for r in json.loads(out_path.read_text())}
    rows, worst = [], 0.0
    for key in sorted(ref):
        a, b = leaves(ref[key]), leaves(got[key])
        diffs = {}
        for k in sorted(set(a) | set(b)):
            if k in ("seconds",) or k.startswith("audit."):
                continue
            x, y = a.get(k), b.get(k)
            if isinstance(x, (int, float)) and isinstance(y, (int, float)) and not isinstance(x, bool):
                d = 0.0 if (isinstance(x, float) and isinstance(y, float)
                            and math.isnan(x) and math.isnan(y)) else abs(float(x) - float(y))
            else:
                d = 0.0 if x == y else float("inf")
            if d != 0.0:
                diffs[k] = (x, y, d)
        m = max((d for _, _, d in diffs.values()), default=0.0)
        worst = max(worst, m)
        rows.append(dict(panel=key[0], fold=key[1], cos_agg=got[key]["cos_agg"], pds=got[key]["pds"],
                         ref_cos_agg=ref[key]["cos_agg"], ref_pds=ref[key]["pds"],
                         n_leaves_compared=len(set(a) | set(b)), max_abs_diff=m,
                         differing=diffs,
                         h1_per_unit_all_rows=got[key]["sizing"]["h1_per_unit_all_rows"],
                         cd4_per_unit_all_rows=got[key]["sizing"]["cd4_per_unit_all_rows"]))
        print(f"[{key[0]}/{key[1]}] cos_agg {got[key]['cos_agg']:.15g} pds {got[key]['pds']:.15g} "
              f"vs 233: max|d| {m:.1e} over {len(set(a) | set(b))} leaves"
              f"{'  DIFFER ' + str(list(diffs)[:5]) if diffs else ''}", flush=True)
    out = dict(rows=rows, worst_abs_diff=worst, n_rows=len(rows),
               sizing_default_named="legacy_all_rows",
               seconds=round(time.time() - t0, 1),
               peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3)
    (WORK / "ka233.json").write_text(json.dumps(out, indent=1, default=str))
    print(f"KA-233: {len(rows)} BASE rows, worst |d| {worst:.1e}, {out['seconds']}s, "
          f"peak {out['peak_rss_gib']:.2f} GiB -> {'PASS' if worst == 0.0 else 'FAIL'}")
    if worst != 0.0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
