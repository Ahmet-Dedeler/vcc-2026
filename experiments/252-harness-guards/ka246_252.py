"""252 known-answer row 1b -- 246's bridge rows through the guarded harness, called EXACTLY as
`experiments/250-deflation-october/bridge250.py` calls it (its `stage246` -> `plolo.build_fold`
on 246's own cache copies, `whiten_with` -> `plolo.Operator(...).apply`, `arms` ->
`plolo.solve_fraction` + 246's deflation/controls + `resid218.score_effects`), part (a) only
(062's held operators / the production operator), one fold per process.  Every row `arms()`
tags with `ka=` is asserted against 246's json at 1e-9 inside `arms()` itself; this driver
records the max |d| per row (expected 0.0e+00) for the five the mandate names:

    STACK   B2_cd4_svd5_resolved   B2_cd4_randE5 (random control)
    B2_cd4_svd5_SHUFDEF (shuffled projector)   C2_partial@0.34

Writes only `work/ka246-<fold>.json`; bridge250's own assertion that nothing lands in 246's
cache directory runs inside `stage246`.

    /usr/bin/time -l uv run python ka246_252.py --fold k562_gwps   (orion_hct116, orion_hek293t, kolf2)
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/250-deflation-october"))

import bridge250 as B  # noqa: E402  (imports lib250 -> bind250/panels233/score233/fold233/plolo)
import lib250 as L  # noqa: E402
import plolo  # noqa: E402

FIVE = ("STACK", "B2_cd4_svd5_resolved", "B2_cd4_randE5", "B2_cd4_svd5_SHUFDEF", "C2_partial@0.34")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True, choices=("k562_gwps", "orion_hct116",
                                                     "orion_hek293t", "kolf2"))
    held = ap.parse_args().fold
    t0 = time.time()
    ref = B.ref246(held)
    base = B.stage246(held)                       # 246's instrument; asserts 246's cache untouched
    part = L.cd4_partial("abc", base["targets"], cd4_root=B.SHIPPED_CD4)
    op246 = B.PROD_OP if held == "kolf2" else B.HELD_OPS / f"held_{held}.npz"
    st = B.whiten_with(base, op246)
    rows = B.arms(st, ref, "op062" if held != "kolf2" else "production", held, part)
    ka = {r["arm"]: r["ka"] for r in rows if "ka" in r}
    five = {a: ka.get(a) for a in FIVE}
    worst = max((v["max_abs_diff"] for v in ka.values()), default=float("nan"))
    out = dict(fold=held, n_targets=len(base["targets"]), operator=str(op246),
               op_name=plolo.registered_operator_name(op246),
               cache=base["fold"]["audit"].get("cache"),
               ka_rows=ka, five=five, worst_abs_diff=worst,
               rows=[{k: v for k, v in r.items() if k != "audit"} for r in rows],
               seconds=round(time.time() - t0, 1),
               peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3)
    (WORK / f"ka246-{held}.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"[{held}] {len(ka)} KA rows vs 246, worst |d| {worst:.1e}; five named: "
          + ", ".join(f"{a}={v['max_abs_diff']:.1e}" if v else f"{a}=ABSENT" for a, v in five.items())
          + f"; {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB")
    if any(v is None or not v["ok"] for v in five.values()) or worst != 0.0:
        print("KA-246 FAIL")
        raise SystemExit(1)
    print("KA-246 PASS")


if __name__ == "__main__":
    main()
