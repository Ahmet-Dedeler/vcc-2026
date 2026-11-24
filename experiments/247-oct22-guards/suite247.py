"""247 — the builder's own regression suite, test by test, under BOTH builders.

`scripts/test_make_submission_pm_cd4.py`'s runner stops at the first failure, and one test
(`test_panel_exclusion_for_operator_and_cd4_residual`) has unpacked two values from a function
that returns three since 222 -- so the suite as a whole has not passed since then.  This runs
every test in isolation under the working-tree builder AND under the pre-247 byte copy (loaded
with `__file__` pinned to `scripts/`, then swapped into the test module as `M`), and reports the
DIFFERENCE: any test that passes under pre-247 and fails under 247 is a regression of this
experiment; anything failing under both predates it.

    uv run python experiments/247-oct22-guards/suite247.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))


def load(path: Path, name: str, pinned: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(pinned)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def run_all(T, label: str) -> dict[str, dict]:
    out = {}
    for test in T.TESTS:
        t0 = time.time()
        try:
            test()
            out[test.__name__] = {"ok": True, "seconds": round(time.time() - t0, 1)}
        except Exception as exc:  # noqa: BLE001
            out[test.__name__] = {"ok": False, "seconds": round(time.time() - t0, 1),
                                  "error": f"{type(exc).__name__}: {exc}"[:300]}
        print(f"[{label}] {'PASS' if out[test.__name__]['ok'] else 'FAIL'} {test.__name__} "
              f"({out[test.__name__]['seconds']}s)"
              + ("" if out[test.__name__]["ok"] else f"  {out[test.__name__]['error'][:120]}"), flush=True)
    return out


def main() -> None:
    T = load(ROOT / "scripts/test_make_submission_pm_cd4.py", "builder_tests_247",
             ROOT / "scripts/test_make_submission_pm_cd4.py")
    new = T.M
    old = load(HERE / "work/pre247/make_submission_pm_cd4.py", "builder_pre247_for_suite",
               ROOT / "scripts/make_submission_pm_cd4.py")
    res_new = run_all(T, "247")
    T.M = old
    res_old = run_all(T, "pre-247")
    T.M = new
    regressions = [k for k in res_new if not res_new[k]["ok"] and res_old.get(k, {}).get("ok")]
    fixed = [k for k in res_new if res_new[k]["ok"] and not res_old.get(k, {}).get("ok")]
    both_fail = [k for k in res_new if not res_new[k]["ok"] and not res_old.get(k, {}).get("ok")]
    summary = {
        "n_tests": len(T.TESTS),
        "pass_247": sum(r["ok"] for r in res_new.values()),
        "pass_pre247": sum(r["ok"] for r in res_old.values()),
        "regressions_247_vs_pre247": regressions,
        "fixed_by_247": fixed,
        "failing_under_both_predate_247": {k: res_new[k]["error"] for k in both_fail},
        "peak_rss_gib": round(new.peak_rss_gib(), 3),
    }
    (WORK / "suite247.json").write_text(json.dumps({"summary": summary, "247": res_new, "pre247": res_old}, indent=1))
    print(json.dumps(summary, indent=1))
    print(f"\nSUITE {'PASS (no regression)' if not regressions else 'FAIL (regression)'}")
    raise SystemExit(1 if regressions else 0)


if __name__ == "__main__":
    main()
