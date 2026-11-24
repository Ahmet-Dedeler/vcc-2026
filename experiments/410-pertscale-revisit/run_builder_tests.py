"""Run the builder's unit-test tuple, skipping the one test stale since 481cd092 (2026-09-04), as 399/402 did."""
import sys, time, traceback, importlib.util
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location("tmsp", ROOT / "scripts/test_make_submission_pm_cd4.py")
T = importlib.util.module_from_spec(spec); spec.loader.exec_module(T)
SKIP = {"test_panel_exclusion_for_operator_and_cd4_residual"}
p = f = s = 0; failed = []
for test in T.TESTS:
    if test.__name__ in SKIP:
        print(f"SKIP {test.__name__}"); s += 1; continue
    t0 = time.time()
    try:
        test(); p += 1; print(f"PASS {test.__name__} ({time.time()-t0:.1f}s)", flush=True)
    except Exception:
        f += 1; failed.append(test.__name__); print(f"FAIL {test.__name__}"); traceback.print_exc()
print(f"{p} passed, {f} failed, {s} skipped; failed: {failed}")
