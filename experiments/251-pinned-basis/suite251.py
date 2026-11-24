"""251 -- the builder's own regression suite, test by test, under BOTH builders (249's
`suite249.py`, byte copy swapped for this experiment's).  Reports regressions of 251 only.

    uv run python experiments/251-pinned-basis/suite251.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
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
    T = load(ROOT / "scripts/test_make_submission_pm_cd4.py", "builder_tests_251",
             ROOT / "scripts/test_make_submission_pm_cd4.py")
    new = T.M
    old = load(WORK / "make_submission_pm_cd4.PRE.py", "builder_pre251_for_suite",
               ROOT / "scripts/make_submission_pm_cd4.py")
    res_new = run_all(T, "251")
    T.M = old
    res_old = run_all(T, "pre-251")
    T.M = new
    regressions = [k for k in res_new if not res_new[k]["ok"] and res_old.get(k, {}).get("ok")]
    fixed = [k for k in res_new if res_new[k]["ok"] and not res_old.get(k, {}).get("ok")]
    both_fail = [k for k in res_new if not res_new[k]["ok"] and not res_old.get(k, {}).get("ok")]
    summary = {
        "n_tests": len(T.TESTS),
        "pass_251": sum(r["ok"] for r in res_new.values()),
        "pass_pre251": sum(r["ok"] for r in res_old.values()),
        "regressions_251_vs_pre251": regressions,
        "fixed_by_251": fixed,
        "failing_under_both_predate_251": {k: res_new[k]["error"] for k in both_fail},
        "peak_rss_gib": round(new.peak_rss_gib(), 3),
    }
    (WORK / "suite251.json").write_text(json.dumps({"summary": summary, "251": res_new, "pre251": res_old}, indent=1))
    print(json.dumps(summary, indent=1))
    print(f"\nSUITE {'PASS (no regression)' if not regressions else 'FAIL (regression)'}")
    raise SystemExit(1 if regressions else 0)


if __name__ == "__main__":
    main()
